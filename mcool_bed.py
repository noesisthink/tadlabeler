#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAD边界检测综合脚本 - 集成DI、SDI和IS三种方法
纯 h5py + numpy + pandas 实现，无需 cooler / pybedtools / bedtools，Windows 可用。
用法: python tad_detection.py -i <mcool文件> -r <分辨率> [选项]
"""

import argparse
import os
import sys
import time

import numpy as np
import h5py
import pandas as pd
from tqdm import tqdm
from hmmlearn import hmm



# ═══════════════════════════════════════════════
# mcool 读取层（替代 cooler）
# ═══════════════════════════════════════════════

import h5py
import pandas as pd
import numpy as np
import time


class McoolReader:
    """高性能 .mcool 读取器：支持染色体筛选与切片读取。"""

    def __init__(self, cool_file: str, resolution: int, chroms: list = None):
        self.cool_file = cool_file
        self.resolution = self._resolve_resolution(resolution)
        self._root_path = f"resolutions/{self.resolution}"

        # 1. 加载基础元数据
        self._bins_df = self._load_bins()

        # 2. 确定目标染色体列表（并排序）
        all_available = self._bins_df["chrom"].unique().tolist()
        if chroms:
            # 只保留文件中存在的染色体，并按自然顺序排序
            self._target_chroms = sorted(
                [c for c in chroms if c in all_available],
                key=lambda c: self._chrom_sort_key(c)
            )
        else:
            self._target_chroms = sorted(
                all_available,
                key=lambda c: self._chrom_sort_key(c)
            )

        print(f" McoolReader 初始化完成 | 分辨率: {self.resolution} | 目标染色体: {len(self._target_chroms)} 条")

    # ── 内部辅助 ──────────────────────────────

    @staticmethod
    def _chrom_sort_key(name: str):
        n = str(name).replace("chr", "")
        if n.isdigit(): return (0, int(n))
        return (1, {"X": 100, "Y": 101, "M": 102, "MT": 102}.get(n.upper(), 200))

    def _resolve_resolution(self, resolution: int) -> int:
        with h5py.File(self.cool_file, "r") as f:
            available = [int(r) for r in f["resolutions"].keys() if r.isdigit()]
        if resolution in available: return resolution
        closest = min(available, key=lambda x: abs(x - resolution))
        return closest

    def _load_bins(self) -> pd.DataFrame:
        """加载 bins 信息并预计算每个染色体的 bin 范围。"""
        with h5py.File(self.cool_file, "r") as f:
            grp = f[self._root_path]
            chrom_ids = grp["bins"]["chrom"][:]
            starts = grp["bins"]["start"][:]
            ends = grp["bins"]["end"][:]
            chrom_names = [c.decode() if isinstance(c, bytes) else c for c in grp["chroms"]["name"][:]]

        chrom_col = [chrom_names[i] for i in chrom_ids]
        return pd.DataFrame({"chrom": chrom_col, "start": starts, "end": ends})

    # ── 公开接口 ──────────────────────────────

    @property
    def chromnames(self) -> list:
        """返回本次任务需要处理的染色体列表"""
        return self._target_chroms

    def bins_for(self, chrom: str) -> pd.DataFrame:
        return self._bins_df[self._bins_df["chrom"] == chrom].reset_index(drop=True)

    def matrix(self, chrom: str) -> np.ndarray:
        """
        利用索引（bin1_offset）实现 O(1) 内存开销的切片读取。
        不再一次性读取整个 pixels 数组。
        """
        print(f"  提取矩阵: {chrom} ...")
        start_time = time.time()

        # 1. 获取该染色体对应的全局 bin ID 范围
        chrom_bins = self._bins_df[self._bins_df["chrom"] == chrom]
        if chrom_bins.empty:
            raise ValueError(f"染色体 {chrom} 不存在")

        id_min = int(chrom_bins.index[0])
        id_max = int(chrom_bins.index[-1])
        n = id_max - id_min + 1

        with h5py.File(self.cool_file, "r") as f:
            root = f[self._root_path]
            pix_grp = root["pixels"]
            # 关键：获取 bin1 的偏移索引，实现精准定位
            bin1_offsets = root["indexes"]["bin1_offset"][:]

            # 2. 确定该染色体在 pixels 数组中的物理起始位置
            idx_start = bin1_offsets[id_min]
            idx_end = bin1_offsets[id_max + 1]

            # 3. 仅读取该染色体相关的像素块（极大节省内存！）
            b1 = pix_grp["bin1_id"][idx_start:idx_end]
            b2 = pix_grp["bin2_id"][idx_start:idx_end]
            cnt = pix_grp["count"][idx_start:idx_end]

        # 4. 构建密集矩阵
        # 注意：由于 .cool 是对称矩阵存储的上三角，我们只需要过滤 b2 越界的情况
        mask = (b2 >= id_min) & (b2 <= id_max)
        row = b1[mask] - id_min
        col = b2[mask] - id_min
        data = cnt[mask].astype(float)

        mat = np.zeros((n, n), dtype=float)
        mat[row, col] = data
        mat[col, row] = data  # 填充对称部分

        print(f"  完成 | 大小: {n}x{n} | 耗时: {time.time() - start_time:.2f}s")
        return mat


# ═══════════════════════════════════════════════
# 公共工具
# ═══════════════════════════════════════════════

def fmt_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


# ═══════════════════════════════════════════════
# 方法一：DI (Directionality Index)
# ═══════════════════════════════════════════════

def calculate_DI(hic_matrix: np.ndarray, window_size: int = 200) -> np.ndarray:
    print("  计算 Directionality Index (DI)...")
    n  = hic_matrix.shape[0]
    DI = np.zeros(n)
    for i in tqdm(range(window_size, n - window_size), desc="  DI"):
        A = np.sum(hic_matrix[i, i - window_size:i])
        B = np.sum(hic_matrix[i, i + 1:i + window_size + 1])
        E = (A + B) / 2
        if E == 0:
            DI[i] = 0
        else:
            sign  = (B - A) / abs(B - A) if B != A else 0
            DI[i] = sign * (((A - E) ** 2 / E) + ((B - E) ** 2 / E))
    return DI


def detect_boundaries_DI(DI_values: np.ndarray, resolution: int) -> list:
    print("  用 HMM 检测 TAD 边界...")
    vals = DI_values.reshape(-1, 1).copy()
    vals = np.nan_to_num(vals)
    vals[np.isinf(vals)] = 0
    model  = hmm.GaussianHMM(n_components=3, covariance_type="full", n_iter=100)
    model.fit(vals)
    states = model.predict(vals)
    bds    = [i * resolution for i in range(1, len(states)) if states[i] != states[i - 1]]
    print(f"  检测到 {len(bds)} 个TAD边界")
    return bds


def run_DI(reader: McoolReader, output_file: str):
    print("\n[1/4] 运行 DI 方法...")
    res    = reader.resolution
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat = reader.matrix(chrom)
            DI  = calculate_DI(mat)
            bds = detect_boundaries_DI(DI, res)
            with open(output_file, "a") as f:
                for b in bds:
                    f.write(f"{fmt_chrom(chrom)}\t{b - res//2}\t{b + res//2}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  DI 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 方法二：SDI (Simplified Directionality Index)
# ═══════════════════════════════════════════════

def calculate_SDI(hic_matrix: np.ndarray,
                  window_sizes: list = None) -> np.ndarray:
    if window_sizes is None:
        window_sizes = [7, 8, 9, 10, 11, 12]
    print("  计算 sDI（多窗口投票）...")
    n      = hic_matrix.shape[0]
    counts = np.zeros(n)
    for w in window_sizes:
        for i in tqdm(range(w, n - w), desc=f"  sDI window={w}"):
            A = np.sum(hic_matrix[i - w:i, i - w:i])
            B = np.sum(hic_matrix[i:i + w, i:i + w])
            if A > B:
                counts[i] -= 1
            elif A < B:
                counts[i] += 1
    return np.where(counts < 0, -1, np.where(counts > 0, 1, 0))


def detect_boundaries_SDI(sdi_values: np.ndarray, resolution: int) -> list:
    bds = [i * resolution for i in range(1, len(sdi_values))
           if sdi_values[i - 1] != sdi_values[i] and sdi_values[i] != 0]
    print(f"  检测到 {len(bds)} 个TAD边界")
    return bds


def run_SDI(reader: McoolReader, output_file: str):
    print("\n[2/4] 运行 SDI 方法...")
    res    = reader.resolution
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat = reader.matrix(chrom)
            sdi = calculate_SDI(mat)
            bds = detect_boundaries_SDI(sdi, res)
            with open(output_file, "a") as f:
                for b in bds:
                    f.write(f"{fmt_chrom(chrom)}\t{b - res//2}\t{b + res//2}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  SDI 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 方法三：IS (Insulation Score)
# ═══════════════════════════════════════════════

def calculate_IS(hic_matrix: np.ndarray, window_size: int) -> np.ndarray:
    print(f"  计算 Insulation Score（窗口={window_size} bins）...")
    n  = hic_matrix.shape[0]
    IS = np.zeros(n)
    for i in tqdm(range(window_size, n - window_size), desc="  IS"):
        IS[i] = np.nansum(
            hic_matrix[i - window_size:i + window_size,
                        i - window_size:i + window_size]
        )
    return IS


def detect_boundaries_IS(is_values: np.ndarray, threshold: float) -> list:
    print(f"  检测 TAD 边界（阈值={threshold}）...")
    std = (is_values - np.nanmean(is_values)) / (np.nanstd(is_values) + 1e-12)
    bds = [i for i in range(1, len(std) - 1)
           if std[i] < std[i - 1] and std[i] < std[i + 1] and std[i] < threshold]
    print(f"  检测到 {len(bds)} 个TAD边界")
    return bds


def run_IS(reader: McoolReader, threshold: float, output_file: str):
    print("\n[3/4] 运行 IS 方法...")
    res         = reader.resolution
    window_bins = 500_000 // res
    print(f"  窗口大小: {window_bins} bins (500kb)")
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat     = reader.matrix(chrom)
            bins    = reader.bins_for(chrom)
            is_vals = calculate_IS(mat, window_bins)
            bds     = detect_boundaries_IS(is_vals, threshold)
            with open(output_file, "a") as f:
                for idx in bds:
                    if idx < len(bins):
                        row = bins.iloc[idx]
                        f.write(f"{fmt_chrom(row['chrom'])}\t{row['start']}\t{row['end']}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  IS 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 快速模式：每条染色体只加载一次矩阵，三方法共用
# ═══════════════════════════════════════════════

def run_all_methods(reader: McoolReader, is_threshold: float,
                    di_file: str, sdi_file: str, is_file: str):
    """
    每条染色体只读取一次矩阵，依次运行 DI / SDI / IS，
    结果分别追加到各自的输出文件。
    原有的 run_DI / run_SDI / run_IS 函数保持不变。
    """
    print("\n[快速模式] 每条染色体单次加载，依次运行 DI / SDI / IS ...")
    res         = reader.resolution
    window_bins = 500_000 // res
    chroms      = reader.chromnames

    # 清空三个输出文件
    for path in (di_file, sdi_file, is_file):
        open(path, "w").close()

    ok = 0
    t_total = time.time()

    for idx, chrom in enumerate(chroms):
        print(f"\n{'─'*50}")
        print(f"  [{idx+1}/{len(chroms)}] 染色体 {chrom}")
        print(f"{'─'*50}")

        try:
            # ── 一次性加载矩阵 ──
            mat  = reader.matrix(chrom)
            bins = reader.bins_for(chrom)

            # ── DI ──
            print("  → DI")
            DI  = calculate_DI(mat)
            bds = detect_boundaries_DI(DI, res)
            with open(di_file, "a") as f:
                for b in bds:
                    f.write(f"{fmt_chrom(chrom)}\t{b - res//2}\t{b + res//2}\n")

            # ── SDI ──
            print("  → SDI")
            sdi = calculate_SDI(mat)
            bds = detect_boundaries_SDI(sdi, res)
            with open(sdi_file, "a") as f:
                for b in bds:
                    f.write(f"{fmt_chrom(chrom)}\t{b - res//2}\t{b + res//2}\n")

            # ── IS ──
            print("  → IS")
            is_vals = calculate_IS(mat, window_bins)
            bds     = detect_boundaries_IS(is_vals, is_threshold)
            with open(is_file, "a") as f:
                for i in bds:
                    if i < len(bins):
                        row = bins.iloc[i]
                        f.write(f"{fmt_chrom(row['chrom'])}\t{row['start']}\t{row['end']}\n")

            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")

    print(f"\n[快速模式] 完成，成功处理 {ok}/{len(chroms)} 条染色体，"
          f"总用时: {time.time()-t_total:.1f}s")
    print(f"  DI  → {di_file}")
    print(f"  SDI → {sdi_file}")
    print(f"  IS  → {is_file}")


# ═══════════════════════════════════════════════
# 步骤四：并集 - 交集（纯 pandas，替代 pybedtools/bedtools）
# ═══════════════════════════════════════════════

def load_bed(path: str) -> pd.DataFrame:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end"], usecols=[0, 1, 2])
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    return df


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """等价于 bedtools merge。"""
    if df.empty:
        return df.copy()
    out = []
    for chrom, grp in df.groupby("chrom", sort=False):
        grp = grp.sort_values("start")
        cs, ce = int(grp.iloc[0]["start"]), int(grp.iloc[0]["end"])
        for _, row in grp.iloc[1:].iterrows():
            if int(row["start"]) <= ce:
                ce = max(ce, int(row["end"]))
            else:
                out.append((chrom, cs, ce))
                cs, ce = int(row["start"]), int(row["end"])
        out.append((chrom, cs, ce))
    return pd.DataFrame(out, columns=["chrom", "start", "end"])


def _has_overlap(as_, ae, bs, be, f: float) -> bool:
    inter = max(0, min(ae, be) - max(as_, bs))
    if inter == 0:
        return False
    return (inter / (ae - as_) >= f) and (inter / (be - bs) >= f)


def intersect_rows(df_a: pd.DataFrame, df_b: pd.DataFrame,
                   f: float) -> pd.DataFrame:
    """返回 df_a 中与 df_b 至少有 f 比例重叠的行（等价 bedtools intersect -u）。"""
    if df_a.empty or df_b.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    b_by_chrom = {c: g.values for c, g in df_b.groupby("chrom")}
    hits = []
    for row in df_a.itertuples(index=False):
        grp = b_by_chrom.get(row.chrom)
        if grp is None:
            continue
        for rb in grp:                      # rb: (chrom, start, end)
            if _has_overlap(row.start, row.end, int(rb[1]), int(rb[2]), f):
                hits.append((row.chrom, row.start, row.end))
                break
    return pd.DataFrame(hits, columns=["chrom", "start", "end"])


def subtract_rows(df_a: pd.DataFrame, df_b: pd.DataFrame,
                  f: float) -> pd.DataFrame:
    """返回 df_a 中 *不* 与 df_b 重叠的行（等价 bedtools intersect -v）。"""
    if df_b.empty:
        return df_a.copy()
    overlap_set = set(
        map(tuple, intersect_rows(df_a, df_b, f)[["chrom", "start", "end"]].values)
    )
    mask = [tuple(r) not in overlap_set
            for r in df_a[["chrom", "start", "end"]].values]
    return df_a[mask].reset_index(drop=True)


def run_intersection(di_file: str, sdi_file: str, is_file: str,
                     output_file: str, overlap: float):
    print("\n[4/4] 计算三种方法的并集减去交集（纯Python）...")

    di_df  = load_bed(di_file)
    sdi_df = load_bed(sdi_file)
    is_df  = load_bed(is_file)

    # 1. 并集 + merge
    union_df = merge_intervals(
        pd.concat([di_df, sdi_df, is_df], ignore_index=True)
    )
    print(f"  并集区域数: {len(union_df)}")

    # 2. 三者交集：DI∩SDI∩IS
    di_sdi_df = intersect_rows(di_df, sdi_df, overlap)
    inter_df  = intersect_rows(di_sdi_df, is_df, overlap)
    print(f"  三者交集区域数: {len(inter_df)}")

    # 3. 并集 - 交集
    result_df = subtract_rows(union_df, inter_df, overlap)
    print(f"  并集减交集区域数: {len(result_df)}")
    result_df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"  结果写入: {output_file}")

    # 4. 各方法独有边界
    di_sdi_ov = intersect_rows(di_df,  sdi_df, overlap)
    di_is_ov  = intersect_rows(di_df,  is_df,  overlap)
    sdi_is_ov = intersect_rows(sdi_df, is_df,  overlap)
    di_only   = subtract_rows(subtract_rows(di_df,  di_sdi_ov, overlap), di_is_ov,  overlap)
    sdi_only  = subtract_rows(subtract_rows(sdi_df, di_sdi_ov, overlap), sdi_is_ov, overlap)
    is_only   = subtract_rows(subtract_rows(is_df,  di_is_ov,  overlap), sdi_is_ov, overlap)

    print("\n===== 统计信息 =====")
    print(f"  DI  边界数: {len(di_df)}")
    print(f"  SDI 边界数: {len(sdi_df)}")
    print(f"  IS  边界数: {len(is_df)}")
    print(f"  三者交集  : {len(inter_df)}")
    print(f"  并集-交集 : {len(result_df)}")
    print(f"  DI  独有  : {len(di_only)}")
    print(f"  SDI 独有  : {len(sdi_only)}")
    print(f"  IS  独有  : {len(is_only)}")
    print("====================")


# ═══════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="TAD边界检测（DI + SDI + IS），纯h5py/pandas实现，Windows可用",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例: python tad_detection.py -i sample.mcool -r 10000 -t 0.1 -o 0.7 -b myhash"
    )
    p.add_argument("-i", "--input",        required=True,  help="输入 .mcool 文件路径")
    p.add_argument("-r", "--resolution",   required=True,  type=int,   help="分析分辨率 (bp)")
    p.add_argument("-t", "--is-threshold", default=0.5,    type=float, help="IS方法阈值 (默认: 0.5)")
    p.add_argument("-o", "--overlap",      default=0.5,    type=float, help="交集重叠度 0-1 (默认: 0.5)")
    p.add_argument("-b", "--hash",         default="",                 help="输出文件名标识符")
    p.add_argument("--fast",  action="store_true",
                   help="快速模式：每条染色体只加载一次矩阵，三方法共用（推荐）")
    return p.parse_args()


def detect_TAD_boundaries(
        mcool_file: str,
        resolution: int,
        selected_chroms: list = None,  # <--- 新增参数：用户选中的染色体列表
        is_threshold: float = 0.5,
        overlap: float = 0.5,
        hash_tag: str = "",
        fast_mode: bool = True,
        output_dir: str = "./tad_results"
):
    """
    支持按染色体筛选的 TAD 边界检测
    """
    os.makedirs(output_dir, exist_ok=True)

    # 如果没有传入，默认跑全基因组（保持向下兼容）
    if not selected_chroms:
        selected_chroms = []

        # 文件命名加入标识，方便区分不同染色体的运行结果
    chrom_suffix = "_all" if not selected_chroms else f"_{len(selected_chroms)}chroms"

    di_out = f"{output_dir}/DI_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    sdi_out = f"{output_dir}/SDI_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    is_out = f"{output_dir}/IS_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    union_out = f"{output_dir}/TAD_final_res{resolution}_{hash_tag}{chrom_suffix}.bed"

    # 修改 McoolReader 的调用，增加染色体过滤
    # 注意：你需要确保你的 McoolReader 类支持这个参数
    reader = McoolReader(mcool_file, resolution, chroms=selected_chroms)

    if fast_mode:
        # run_all_methods 内部遍历 reader 时将只看到选中的染色体
        run_all_methods(reader, is_threshold, di_out, sdi_out, is_out)
    else:
        run_DI(reader, di_out)
        run_SDI(reader, sdi_out)
        run_IS(reader, is_threshold, is_out)

    # 这里的交集计算逻辑通常是基于 BED 文件的，只要前面的 BED 只有特定染色体，
    # 这里的 union_out 自然也就只有特定染色体
    run_intersection(di_out, sdi_out, is_out, union_out, overlap)

    print(f"\n===== TAD检测完成 ({'Partial' if selected_chroms else 'Full'}) =====")
    return {
        "DI": di_out,
        "SDI": sdi_out,
        "IS": is_out,
        "TAD_final": union_out
    }

def main():
    results = detect_TAD_boundaries(
        mcool_file="/extra_mode/GSE63525_GM12878_insitu_DpnII_combined.mcool",
        resolution=10000,
        is_threshold=0.3,
        overlap=0.7,
        hash_tag="sample1",
        fast_mode=True,
        output_dir="./"
    )

    print(results)


if __name__ == "__main__":
    main()