import os
import winbbi


import pandas as pd
import numpy as np
import h5py
import asyncio
import uuid
import json
from aiohttp import web
import aiohttp_cors
from collections import OrderedDict
from mcool_bed import detect_TAD_boundaries
from sql import init_db, add_token_record, update_tad_outputs, mark_running, mark_failed, get_record
from dactor import test
import io
import sys

# 核心修复：强制标准输出和错误流使用 UTF-8 编码
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# =========================================================
# 1. BED 渲染引擎 (智能解析版)
# =========================================================
class BedTileRenderer:
    def __init__(self, bed_path, tile_size=1024,prediction_path=None):
        self.path = bed_path
        self.tile_size = tile_size
        self.prediction_map = {}
        try:
            # 1. 智能读取：sep=r'\s+' 兼容空格和 Tab，comment='#' 跳过注释
            # 先不指定 dtype，防止表头导致 int 转换崩溃
            self.df = pd.read_csv(
                bed_path,
                sep=r'\s+',
                header=None,
                comment='#',
                engine='python'
            )

            # 2. 自动识别列数并处理
            col_count = len(self.df.columns)
            if col_count >= 4:
                self.df = self.df.iloc[:, :4]
                self.df.columns = ['chrom', 'start', 'end', 'name']
            elif col_count == 3:
                self.df = self.df.iloc[:, :3]
                self.df.columns = ['chrom', 'start', 'end']
                self.df['name'] = ""  # 补齐列
            else:
                raise ValueError(f"BED 文件至少需要 3 列 (chrom, start, end)，检测到 {col_count} 列")

            # 3. 强制类型清洗 (关键：处理掉可能存在的 'start', 'end' 字符串表头)
            self.df['start'] = pd.to_numeric(self.df['start'], errors='coerce')
            self.df['end'] = pd.to_numeric(self.df['end'], errors='coerce')

            # 删除无法转换为数字的行（即表头行或脏数据）
            self.df = self.df.dropna(subset=['start', 'end'])

            self.df['start'] = self.df['start'].astype(int)
            self.df['end'] = self.df['end'].astype(int)
            self.df['chrom'] = self.df['chrom'].astype(str).str.strip().str.replace('^chr', '', case=False, regex=True)
            self.df = self.df.fillna("")

            # 4. 建立索引映射
            self.chrom_groups = {name: group for name, group in self.df.groupby('chrom')}

            # 这里的分辨率需要与你的 Hi-C mcool 保持一致或自定义层级
            self.resolutions = [1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]

            print(f" BED 加载成功: {bed_path}, 解析列数: {col_count}, 总记录: {len(self.df)}")
            if prediction_path and os.path.exists(prediction_path):
                pred_df = pd.read_csv(prediction_path, sep=r'\s+', header=0)
                pred_df['chrom'] = pred_df['chrom'].astype(str).str.replace('^chr', '', regex=True)
                for _, row in pred_df.iterrows():
                    key = (str(row['chrom']), int(row['start']), int(row['end']))
                    self.prediction_map[key] = float(row['prediction'])
                print(f" 预测概率加载成功: {len(self.prediction_map)} 条")
        except Exception as e:
            print(f" BED 加载失败: {str(e)}")
            raise e

    def _normalize_chrom(self, chrom):
        return str(chrom).strip().replace('chr', '').replace('Chr', '')

    def get_all_chroms(self):
        chrom_info = []
        offset = 0
        # 按照 DataFrame 中的染色体顺序获取
        for name in self.df['chrom'].unique():
            group = self.chrom_groups[name]
            length = int(group['end'].max())
            chrom_info.append({"name": name, "length": length, "offset": offset})
            offset += length
        return chrom_info

    def get_info(self, chrom):
        chrom_key = self._normalize_chrom(chrom)
        if chrom_key not in self.chrom_groups:
            return {"error": f"未找到染色体: {chrom_key}", "available": list(self.chrom_groups.keys())}
        group = self.chrom_groups[chrom_key]
        return {
            "chromosome": chrom_key,
            "max_zoom": len(self.resolutions) - 1,
            "available_resolutions": self.resolutions,
            "genome_length": int(group['end'].max()),
            "tile_size": self.tile_size
        }

    async def fetch_tile(self, chrom, zoom, x):
        chrom_key = self._normalize_chrom(chrom)
        if chrom_key not in self.chrom_groups: return []

        # 确保 zoom 不越界
        z = min(max(0, zoom), len(self.resolutions) - 1)
        res = self.resolutions[z]

        start_bp, end_bp = x * self.tile_size * res, (x + 1) * self.tile_size * res
        df = self.chrom_groups[chrom_key]

        # 范围查询优化
        subset = df[(df['end'] >= start_bp) & (df['start'] <= end_bp)].copy()

        # 抽样逻辑：防止低倍率下数据量过大导致前端卡死
        if zoom < 3 and len(subset) > 1000:
            subset = subset.sample(n=1000)

        return subset.to_dict(orient='records')




class BigWigRenderer:
    def __init__(self, bw_path):
        self.path = bw_path
        try:
            self.reader = winbbi.BbiFile(bw_path)
            self.header = self.reader.header  # {chrom: length}
            # 保持与 Hi-C 和 BED 一致的分辨率层级
            self.resolutions = [1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]
            print(f" BigWig 加载成功: {bw_path}")
        except Exception as e:
            print(f" BigWig 加载失败: {str(e)}")
            raise e

    def _normalize_chrom(self, chrom):
        return str(chrom).strip().replace('chr', '').replace('Chr', '')

    def get_all_chroms(self):
        """生成全基因组偏移表，用于绝对坐标对齐"""
        chrom_info = []
        offset = 0
        for name, length in self.header.items():
            clean_name = self._normalize_chrom(name)
            chrom_info.append({"name": clean_name, "length": int(length), "offset": offset})
            offset += int(length)
        return chrom_info

    async def fetch_abs_range(self, abs_start, abs_end, res, tile_size=256):
        """
        核心逻辑：输入全局绝对坐标，返回该范围内的信号数组
        """
        # 获取偏移表
        chrom_list = self.get_all_chroms()
        # 初始化结果数组（对应瓦片的每个像素）
        # 注意：BigWig 的 fetch 通常直接返回 bins 个数值
        combined_data = np.zeros(tile_size)

        for item in chrom_list:
            name, length, offset = item['name'], item['length'], item['offset']
            chrom_end_abs = offset + length

            # 检查当前染色体是否在请求的绝对范围内
            if chrom_end_abs >= abs_start and offset <= abs_end:
                # 计算在该染色体内的相对查询范围
                rel_start = max(0, int(abs_start - offset))
                rel_end = min(length, int(abs_end - offset))

                if rel_start >= rel_end: continue

                # 映射到结果数组中的索引位置（保持对齐）
                # 这里我们利用 winbbi 的 fetch(bins=...)
                # 简化处理：对于跨染色体瓦片，我们分段请求
                target_chrom = name if name in self.header else f"chr{name}"

                try:
                    # 计算该片段占瓦片的比例
                    chunk_bins = int(((rel_end - rel_start) / (abs_end - abs_start)) * tile_size)
                    if chunk_bins > 0:
                        vals = self.reader.fetch(target_chrom, rel_start, rel_end, bins=chunk_bins)
                        # 将结果填入对应位置（这里需要精细计算像素偏移，略作简化）
                        # 实际项目中通常 BigWig 不会跨染色体显示，或者前端处理拼接
                        return [float(v) if not np.isnan(v) else 0.0 for v in vals]
                except:
                    pass

        return combined_data.tolist()

    def close(self):
        self.reader.close()

# =========================================================
# 2. HiC 渲染引擎
# =========================================================
class HiCRangeRenderer:
    def __init__(self, mcool_path):
        self.path = mcool_path
        self.f = h5py.File(mcool_path, "r")
        raw_res = [int(k) for k in self.f["resolutions"].keys() if k.isdigit()]
        self.resolutions = sorted(raw_res, reverse=True)
        self.data_map = self._load_metadata()

    def _load_metadata(self):
        data = {}
        for res in self.resolutions:
            g = self.f[f"resolutions/{res}"]
            chrom_names = g["chroms"]["name"][:].astype(str)
            chrom_lengths = g["chroms"]["length"][:]
            chrom_offsets = g["indexes"]["chrom_offset"][:]
            chrom_info = {
                name.replace("chr", ""): {
                    "length": int(length),
                    "start_bin": int(chrom_offsets[i]),
                    "n_bins": int(np.ceil(length / res))
                } for i, (name, length) in enumerate(zip(chrom_names, chrom_lengths))
            }
            data[res] = {
                "chrom_info": chrom_info,
                "pixels": g["pixels"],
                "bin_offsets": g["indexes"]["bin1_offset"][:],
                "weights": g["bins"]["weight"][:] if "weight" in g["bins"] else None
            }
        return data

    def get_all_chroms(self):
        res = self.resolutions[-1]
        info_map = self.data_map[res]["chrom_info"]
        chrom_info = []
        offset = 0
        for name, info in info_map.items():
            chrom_info.append({"name": name, "length": info["length"], "offset": offset})
            offset += info["length"]
        return chrom_info

    async def fetch_range_data(self, res, b1_range_abs, b2_range_abs):
        d = self.data_map.get(res)
        if not d: return [], res

        p = d["pixels"]
        offsets = d["bin_offsets"]
        w = d["weights"]

        total_bins = len(offsets) - 1
        abs_s1, abs_e1 = max(0, min(total_bins, b1_range_abs[0])), max(0, min(total_bins, b1_range_abs[1]))
        abs_s2, abs_e2 = max(0, min(total_bins, b2_range_abs[0])), max(0, min(total_bins, b2_range_abs[1]))

        if abs_s1 >= abs_e1: return [], res

        idx_start, idx_end = int(offsets[abs_s1]), int(offsets[abs_e1])
        if idx_start >= idx_end: return [], res

        bin2_ids_abs = p["bin2_id"][idx_start:idx_end]
        counts = p["count"][idx_start:idx_end]

        diffs = np.diff(offsets[abs_s1:abs_e1 + 1])
        rows_abs = np.repeat(np.arange(abs_s1, abs_e1), diffs)

        mask = (bin2_ids_abs >= abs_s2) & (bin2_ids_abs < abs_e2)
        if not np.any(mask): return [], res

        f_rows, f_cols, f_vals = rows_abs[mask], bin2_ids_abs[mask], counts[mask].astype(float)

        if w is not None:
            valid = (f_rows < len(w)) & (f_cols < len(w))
            f_rows, f_cols, f_vals = f_rows[valid], f_cols[valid], f_vals[valid]
            norm_vals = w[f_rows.astype(int)] * w[f_cols.astype(int)]
            f_vals = np.nan_to_num(f_vals * norm_vals)

        results = [{"x": int(c * res), "y": int(r * res), "v": float(v)}
                   for r, c, v in zip(f_rows, f_cols, f_vals) if v > 0]
        return results, res

    def close(self):
        self.f.close()

    def get_genome_index(self):
        # 1. 选一个分辨率（通常用最小分辨率获取 chrom 信息）
        res = self.resolutions[-1]
        d = self.data_map[res]

        chrom_info = []
        current_offset = 0

        # 2. 遍历 .mcool 内部存储的染色体顺序
        for name, info in d["chrom_info"].items():
            chrom_info.append({
                "name": name,
                "length": info["length"],
                "offset": current_offset
            })
            current_offset += info["length"]

        return {
            "total_length": current_offset,
            "chromosomes": chrom_info
        }


# =========================================================
# 3. 统一文件管理器
# =========================================================
class GlobalManager:
    def __init__(self):
        self.renderers = {}

    def register(self, path, file_type,prediction_path=None):
        for token, r in self.renderers.items():
            if r.path == path: return token
        token = str(uuid.uuid4())
        self.renderers[token] = BedTileRenderer(path,prediction_path=prediction_path) if file_type == "bed" else HiCRangeRenderer(path)
        return token

    def get(self, token):
        return self.renderers.get(token)

    def close_all(self):
        for r in self.renderers.values():
            if hasattr(r, 'close'): r.close()

    def get_abs_info(self, token):
        renderer = self.get(token)
        if not renderer: return {}, 0, []
        # Hi-C 渲染器逻辑 (h5py)
        if hasattr(renderer, 'f'):
            res = renderer.resolutions[-1]
            chroms = renderer.f[f"resolutions/{res}/chroms"]
            names, lengths = chroms["name"][:].astype(str), chroms["length"][:]
            offsets, current_total, chrom_list = {}, 0, []
            for n, l in zip(names, lengths):
                clean = n.replace("chr", "").replace("Chr", "")
                offsets[clean] = current_total
                chrom_list.append({"name": clean, "length": int(l), "offset": current_total})
                current_total += int(l)
            return offsets, current_total, chrom_list

        # BigWig 渲染器逻辑 (winbbi)
        elif isinstance(renderer, BigWigRenderer):
            chrom_list = renderer.get_all_chroms()
            offsets = {item['name']: item['offset'] for item in chrom_list}
            total_len = sum(item['length'] for item in chrom_list)
            return offsets, total_len, chrom_list

        # BED 渲染器逻辑 (pandas)
        else:
            chrom_list = renderer.get_all_chroms()
            offsets = {item['name']: item['offset'] for item in chrom_list}
            total_len = sum(item['length'] for item in chrom_list)
            return offsets, total_len, chrom_list

    def abs_to_rel(self, abs_pos, offsets):
        """将全局绝对位置转回 (染色体, 相对位置)"""
        last_chrom = None
        for name, offset in offsets.items():
            if abs_pos < offset:
                break
            last_chrom = name
        if last_chrom:
            return last_chrom, abs_pos - offsets[last_chrom]
        return None, 0


# =========================================================
# 4. API 路由处理函数
# =========================================================

async def handle_register(request):
    try:
        # 1. 强制使用 UTF-8 解码，防止 Windows 路径编码问题
        body = await request.read()
        data = json.loads(body.decode('utf-8'))

        path = data.get("path")
        if not path:
            return web.json_response({"error": "path is required"}, status=400)

        # 2. 强力清洗路径（后端双保险）
        path = path.replace('\u202a', '').replace('\u202b', '').strip()

        f_type = data.get("type", "hic")
        prediction_path = data.get("prediction_path", None)
        token = request.app['manager'].register(path, f_type, prediction_path=prediction_path)

        if f_type == "hic":
            # 3. 使用 try-except 包裹数据库操作，防止因为参数类型不匹配导致 500
            try:
                add_token_record(
                    token=token,
                    mcool_path=path,
                    resolution=data.get("resolution"),
                    is_threshold=data.get("is_threshold"),
                    overlap=data.get("overlap"),
                    hash_tag=data.get("hash_tag"),
                    fast_mode=data.get("fast_mode")
                )
            except Exception as db_err:
                print(f"Database error (ignored): {db_err}")
        return web.json_response({"token": token, "type": f_type})
    except json.JSONDecodeError as e:
        return web.json_response({"error": f"Invalid JSON format: {str(e)}"}, status=400)
    except Exception as e:
        return web.json_response({"error": str(e)}, status=500)



async def get_genome_index(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    data = r.get_all_chroms()
    return web.json_response({"total_length": sum(i['length'] for i in data), "chromosomes": data})






# --- BED 路由 ---
async def bed_info(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r:return web.json_response({"error": "Invalid token"}, status=404)
    return web.json_response(r.get_info(req.match_info['chrom']))


async def handle_bed_tile(req):
    token = req.match_info['token']
    # 1. 安全获取参数，防止索引越界
    try:
        z = int(req.match_info['z'])
        x = int(req.match_info['x'])
    except ValueError:
        return web.json_response({"error": "Invalid coordinates"}, status=400)

    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 2. 获取基因组偏移量
    offsets, total_len, _ = req.app['manager'].get_abs_info(token)

    # 3. 边界检查：防止 z 超过 resolutions 数组长度导致 500
    safe_z = min(max(0, z), len(r.resolutions) - 1)
    res = r.resolutions[safe_z]
    tile_size = 256

    # 4. 计算该瓦片的全局范围
    abs_start = x * tile_size * res
    abs_end = (x + 1) * tile_size * res

    results = []

    # 5. 遍历染色体偏移表
    for mcool_chrom_name, offset in offsets.items():
        # --- 核心修复：归一化染色体名称，确保 "chr7" 和 "7" 能匹配 ---
        clean_name = str(mcool_chrom_name).replace('chr', '').replace('Chr', '')

        # 检查 BED 中是否有这条染色体的数据
        df = r.chrom_groups.get(clean_name)
        if df is None or df.empty:
            continue

        # 获取该染色体的最大长度（基于 BED 数据）
        chrom_len = int(df['end'].max())
        chrom_end_abs = offset + chrom_len

        # 6. 判断染色体是否在当前瓦片的绝对范围内
        if chrom_end_abs >= abs_start and offset <= abs_end:

            # 计算在该染色体内的相对查询范围
            rel_start = max(0, abs_start - offset)
            rel_end = abs_end - offset

            # 筛选重叠区间
            subset = df[(df['end'] >= rel_start) & (df['start'] <= rel_end)].copy()

            if not subset.empty:
                # 转换为绝对坐标供前端对齐 Hi-C 矩阵
                # handle_bed_tile 里 subset 转 dict 前加两列
                subset['raw_chrom'] = clean_name
                subset['raw_start'] = subset['start']  # 加偏移前
                subset['raw_end'] = subset['end']
                subset['start'] = subset['start'] + offset
                subset['end'] = subset['end'] + offset

                # 将该段数据加入结果集
                results.extend(subset.to_dict(orient='records'))

    # 7. 返回结果
    return web.json_response({
        "tile_id": f"{z}.{x}",
        "abs_range": [abs_start, abs_end],
        "resolution": res,
        "data": results
    })


async def handle_1d_tile(req):
    token = req.match_info['token']
    try:
        z = int(req.match_info['z'])
        x = int(req.match_info['x'])
    except ValueError:
        return web.json_response({"error": "Invalid coordinates"}, status=400)

    manager = req.app['manager']
    r = manager.get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 基础参数计算
    safe_z = min(max(0, z), len(r.resolutions) - 1)
    res = r.resolutions[safe_z]
    tile_size = 256  # 对应前端每个瓦片的像素宽度
    abs_start = x * tile_size * res
    abs_end = (x + 1) * tile_size * res

    # 获取全局偏移量（用于跨文件对齐）
    offsets, total_len, _ = manager.get_abs_info(token)

    # --- 分支 A: 处理 BigWig (连续信号) ---
    if isinstance(r, BigWigRenderer):
        # 根据绝对坐标反推当前的染色体和相对位置
        chrom_name, rel_pos = manager.abs_to_rel(abs_start, offsets)
        if not chrom_name:
            return web.json_response({"tile_id": f"{z}.{x}", "data": [0.0] * tile_size, "type": "bigwig"})

        # winbbi 兼容性处理
        target_chrom = chrom_name if chrom_name in r.header else f"chr{chrom_name}"
        try:
            # 直接抓取该范围内的 bins
            raw_values = r.reader.fetch(
                target_chrom,
                int(abs_start - offsets[chrom_name]),
                int(abs_end - offsets[chrom_name]),
                bins=tile_size
            )
            data = [float(v) if not np.isnan(v) else 0.0 for v in raw_values]
        except:
            data = [0.0] * tile_size

        return web.json_response({
            "tile_id": f"{z}.{x}",
            "type": "bigwig",
            "data": data,
            "abs_range": [abs_start, abs_end]
        })

    # --- 分支 B: 处理 BED (区间标注) ---
    else:
        results = []
        for mcool_chrom_name, offset in offsets.items():
            clean_name = str(mcool_chrom_name).replace('chr', '').replace('Chr', '')
            df = r.chrom_groups.get(clean_name)
            if df is None or df.empty: continue

            chrom_len = int(df['end'].max())
            chrom_end_abs = offset + chrom_len

            if chrom_end_abs >= abs_start and offset <= abs_end:
                rel_start = max(0, abs_start - offset)
                rel_end = abs_end - offset
                subset = df[(df['end'] >= rel_start) & (df['start'] <= rel_end)].copy()
                if not subset.empty:
                    # 转换为全局绝对坐标返回
                    subset['start'] += offset
                    subset['end'] += offset
                    results.extend(subset.to_dict(orient='records'))

        return web.json_response({
            "tile_id": f"{z}.{x}",
            "type": "bed",
            "data": results,
            "abs_range": [abs_start, abs_end]
        })

async def bed_tiles(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    chrom, zoom, x = req.match_info['chrom'], int(req.match_info['zoom']), int(req.match_info['x'])
    data = await r.fetch_tile(chrom, zoom, x)
    return web.json_response({"data": data, "zoom": zoom, "x": x})


async def handle_mcool_chroms(request):
    """获取 mcool 文件中包含的染色体列表"""
    token = request.match_info['token']
    renderer = request.app['manager'].get(token)

    if not renderer or not hasattr(renderer, 'data_map'):
        return web.json_response({"error": "Invalid token or not a HiC file"}, status=404)

    # 从最低分辨率中获取染色体信息（通常所有分辨率的染色体是一致的）
    res = renderer.resolutions[-1]
    chrom_info = renderer.data_map[res]["chrom_info"]

    # 返回染色体名称列表，例如 ["1", "2", "X", ...]
    return web.json_response({
        "chromosomes": list(chrom_info.keys())
    })

# --- HiC 路由 ---
async def hic_info(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    chrom = req.match_info['chrom'].replace("chr", "")
    res_key = r.resolutions[-1]
    info = r.data_map[res_key]["chrom_info"].get(chrom, {})
    return web.json_response({"length": info.get("length", 0), "available_resolutions": r.resolutions})


async def hic_range(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)

    chrom = req.match_info['chrom'].replace("chr", "")
    res = int(req.query.get("res", r.resolutions[0]))
    s1, e1 = int(req.query.get("start", 0)), int(req.query.get("end", 0))
    s2, e2 = int(req.query.get("start2", s1)), int(req.query.get("end2", e1))

    d = r.data_map.get(res)
    if not d or chrom not in d["chrom_info"]:
        return web.json_response({"error": "Chrom not found"}, status=404)

    offset = d["chrom_info"][chrom]["start_bin"]
    data, actual_res = await r.fetch_range_data(
        res,
        (offset + s1 // res, offset + e1 // res),
        (offset + s2 // res, offset + e2 // res)
    )
    return web.json_response({"resolution": actual_res, "data": data})


async def hic_global_range(req):
    r = req.app['manager'].get(req.match_info.get('token'))
    if not r: return web.json_response({"error": "Invalid token"}, status=404)

    res = int(req.query.get("res", r.resolutions[0]))
    s1, e1 = int(req.query.get("start", 0)), int(req.query.get("end", 0))
    s2, e2 = int(req.query.get("start2", s1)), int(req.query.get("end2", e1))

    data, actual_res = await r.fetch_range_data(res, (s1 // res, e1 // res), (s2 // res, e2 // res))
    return web.json_response({"resolution": actual_res, "data": data})


async def handle_hic_tile(req):
    r = req.app['manager'].get(req.match_info['token'])
    z, x, y = int(req.match_info['z']), int(req.match_info['x']), int(req.match_info['y'])

    tile_size = 256
    res = r.resolutions[z]

    # Hi-C 内部本身就支持绝对 bin 索引查询
    # 我们直接传入全局 bin 范围
    data, actual_res = await r.fetch_range_data(
        res,
        (x * tile_size, (x + 1) * tile_size),
        (y * tile_size, (y + 1) * tile_size)
    )

    return web.json_response({
        "tile_id": f"{z}.{x}.{y}",
        "data": data
    })



async def bed_get_prediction(req):
    """
    GET /api/bed/prediction/{token}?chrom=7&start=65000&end=75000
    返回该区间的预测概率
    """
    token = req.match_info['token']
    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    chrom = req.query.get('chrom', '').replace('chr', '')
    try:
        start = int(req.query.get('start', 0))
        end   = int(req.query.get('end', 0))
    except ValueError:
        return web.json_response({"error": "Invalid coordinates"}, status=400)

    prob = r.prediction_map.get((chrom, start, end), None)
    return web.json_response({
        "chrom": chrom,
        "start": start,
        "end":   end,
        "prediction": prob  # None 表示没有预测值
    })


async def handle_tad_run(request):
    """
    核心修改：只对最终重合的 TAD_final 文件进行预测
    """
    data = await request.json()
    token = data.get('token')
    if not token:
        return web.json_response({"error": "Token is required"}, status=400)

    renderer = request.app['manager'].get(token)
    if not renderer or not hasattr(renderer, 'path'):
        return web.json_response({"error": "Invalid token or not a HiC file"}, status=400)

    # 标记任务开始
    mark_running(token)

    # 获取参数
    record = get_record(token) or {}
    resolution = data.get('resolution', record.get('resolution', 10000))
    is_threshold = data.get('is_threshold', record.get('is_threshold', 0.3))
    overlap = data.get('overlap', record.get('overlap', 0.7))
    hash_tag = data.get('hash_tag', record.get('hash_tag', token[:8]))
    fast_mode = data.get('fast_mode', record.get('fast_mode', True))
    output_dir = data.get('output_dir', './tad_results')

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    selected_chroms = data.get('selected_chroms', [])  # 例如 ["1", "7"]
    try:
        # 1. 运行 TAD 检测，生成包括 TAD_final 在内的所有文件
        output_paths = await asyncio.get_event_loop().run_in_executor(
            None,
            lambda: detect_TAD_boundaries(
                mcool_file=renderer.path,
                selected_chroms=selected_chroms,  # <-- 传递给算法
                resolution=resolution,
                is_threshold=is_threshold,
                overlap=overlap,
                hash_tag=hash_tag,
                fast_mode=fast_mode,
                output_dir=output_dir
            )
        )

        # 2. 只取最终重合的 TAD_final 文件进行预测（核心修改点）
        tad_final_path = output_paths.get('TAD_final')
        prediction_output = None

        if tad_final_path and os.path.exists(tad_final_path):
            print(f" 开始预测最终重合的 TAD 区间：{tad_final_path}")

            # 构造预测输出路径
            prediction_filename = f"prediction_TAD_final_{hash_tag}.txt"
            prediction_path = os.path.join(output_dir, prediction_filename)

            # 运行预测（只预测最终重合的 TAD 区间）
            prediction_output = await asyncio.get_event_loop().run_in_executor(
                None,
                lambda: test.predict_tad_boundary(
                    bed_file_path=tad_final_path,
                    fasta_file_path=data.get('fasta_path', './extra_mode/hg19.fa'),
                    model_path=data.get('model_path', './extra_mode/GM12878_tad_boundary_model.h5'),
                    output_file_path=prediction_path
                )
            )
            print(f" 最终 TAD 预测完成：{prediction_output}")
        else:
            print(f"  未生成 TAD_final 文件，跳过预测")
        # 3. 更新数据库（包含最终预测路径）
        update_tad_outputs(
            token=token,
            output_paths=output_paths,
            prediction_output=prediction_output  # 只存储最终 TAD 的预测结果
        )
        return web.json_response({
            "status": "done",
            "output_paths": output_paths,
            "prediction_output": prediction_output,
            "message": "TAD detection completed (only final overlapping TADs predicted)"
        })

    except Exception as e:
        error_msg = str(e)
        mark_failed(token, error_msg)
        return web.json_response({
            "status": "failed",
            "error": error_msg
        }, status=500)

async def handle_tad_status(request):
    """
    GET /api/tad/status/{token}
    查询 TAD 任务状态及所有入库字段。
    """
    token  = request.match_info['token']
    record = get_record(token)
    if not record:
        return web.json_response({"error": "Not found"}, status=404)
    return web.json_response(record)


# =========================================================
# 5. 服务启动
# =========================================================
async def on_shutdown(app):
    app['manager'].close_all()


# --- 在 handle_functions 区域添加 ---
async def handle_health(request):
    """用于 Electron 探测后端是否启动完成"""
    return web.json_response({"status": "ok", "service": "genomics-api"})

# --- 查询所有 token ---
async def handle_get_all_records(request):
    """
    GET /api/records
    返回数据库中所有 token 记录。
    """
    from sql import get_all_records
    records = get_all_records()
    return web.json_response({"records": records})


async def handle_genome_index(request):
    token = request.match_info['token']
    manager = request.app['manager']
    renderer = manager.get(token)

    if not renderer:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 调用上面的解析逻辑
    index_data = renderer.get_genome_index()
    return web.json_response(index_data)




# --- 删除单条记录 ---
async def handle_delete_record(request):
    """
    DELETE /api/records/{token}
    删除指定 token 的记录。
    """
    from sql import delete_record

    token = request.match_info['token']
    success = delete_record(token)
    if success:
        return web.json_response({"status": "deleted", "token": token})
    else:
        return web.json_response({"error": "Token not found"}, status=404)


async def bed_global_info(req):
    """
    GET /api/bed/global-info/{token}
    返回全局 BED 信息，不需要染色体参数
    供 TAD BED（绝对坐标）使用
    """
    token = req.match_info['token']
    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 获取所有染色体及总长度
    _, total_len, chrom_list = req.app['manager'].get_abs_info(token)

    return web.json_response({
        "genome_length": total_len,
        "available_resolutions": r.resolutions,
        "chromosomes": chrom_list
    })


# --- 删除所有记录 ---
async def handle_delete_all_records(request):
    """
    DELETE /api/records
    删除数据库中所有 token 记录。
    """
    from sql import delete_all_records

    count = delete_all_records()
    return web.json_response({"status": "all_deleted", "deleted_count": count})


def main():
    init_db()   # 建表（幂等，已存在不影响）

    app = web.Application(client_max_size=1024 ** 2 * 100)
    app['manager'] = GlobalManager()
    app.on_shutdown.append(on_shutdown)

    app.add_routes([
        web.get("/api/health", handle_health),  # <--- 添加这一行
        web.post("/api/register",              handle_register),
        web.get('/api/genome/index/{token}',   get_genome_index),
        web.get("/api/records", handle_get_all_records),
        # BED
        web.get('/api/bed/info/{token}/{chrom}',    bed_info),
        web.get('/api/bed/tiles/{token}/{chrom}/{zoom}/{x}', bed_tiles),
        web.delete("/api/records/{token}", handle_delete_record),
        web.delete("/api/records", handle_delete_all_records),
        web.get('/api/v1/tiles/1d/{token}/{z}/{x}', handle_1d_tile),
        web.get('/api/bed/global-info/{token}', bed_global_info),
        web.get('/api/bed/prediction/{token}', bed_get_prediction),
        # HiC
        web.get("/api/mcool/chroms/{token}", handle_mcool_chroms),
        web.get("/api/hic/info/{token}/{chrom}",  hic_info),
        web.get("/api/hic/range/{token}/{chrom}", hic_range),
        web.get("/api/hic/global-range/{token}",  hic_global_range),
        web.get('/api/v1/tiles/2d/{token}/{z}/{x}/{y}', handle_hic_tile),
        # TAD
        web.post("/api/tad/run",               handle_tad_run),
        web.get("/api/tad/status/{token}",     handle_tad_status),
        web.get('/api/genome/index/{token}', handle_genome_index)
    ])
    cors = aiohttp_cors.setup(app, defaults={
        "*": aiohttp_cors.ResourceOptions(allow_headers="*", allow_methods="*")
    })
    for route in list(app.router.routes()):
        cors.add(route)
    print(" 基因组综合数据 API 服务已启动 [Port: 5001]")
    web.run_app(app, port=5001)

if __name__ == "__main__":
    main()