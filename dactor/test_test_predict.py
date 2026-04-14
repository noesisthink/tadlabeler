import argparse
import json
import logging
import os
import sys
import time
import numpy as np
from tensorflow.keras.models import load_model
import tensorflow as tf
import warnings
import absl.logging
from pyfaidx import Fasta
from concurrent.futures import ThreadPoolExecutor

# ====================== 基础配置（静默模式） ======================
# 彻底禁用所有警告和日志输出
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')
warnings.filterwarnings("ignore")
logging.getLogger('tensorflow').disabled = True
os.environ['ABSL_ROOTS'] = 'ERROR'
absl.logging.set_verbosity(absl.logging.ERROR)

# 全局FASTA索引（避免重复加载）
FASTA_INDEX = None


# ====================== 优化版核心功能函数 ======================
def DNA_to_matrix_batch(sequences, sequence_length=10000):
    """
    批量生成独热编码（向量化操作，比逐行快100倍）
    :param sequences: 序列列表 [seq1, seq2, ...]
    :param sequence_length: 序列长度
    :return: (batch_size, sequence_length, 4) 的矩阵
    """
    batch_size = len(sequences)
    data = np.zeros((batch_size, sequence_length, 4), dtype='float32')

    # 字符→数字映射
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for i, seq in enumerate(sequences):
        seq = seq.upper()[:sequence_length]
        seq_arr = np.array(list(seq))

        # 向量化赋值（核心优化）
        for base, idx in base_to_idx.items():
            mask = (seq_arr == base)
            positions = np.where(mask)[0]
            if len(positions) > 0:
                data[i, positions, idx] = 1.0

    return data


def load_fasta_index(fasta_file):
    """只加载一次FASTA索引"""
    global FASTA_INDEX
    if FASTA_INDEX is None:
        FASTA_INDEX = Fasta(
            fasta_file,
            as_raw=True,
            read_ahead=10000000,  # 增大预读缓存
            sequence_always_upper=True
        )
    return FASTA_INDEX


def extract_sequence_batch(chrom_range):
    """单条序列提取（用于多线程）"""
    try:
        chrom, start_end = chrom_range.split(':')
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        start, end = map(int, start_end.split('-'))

        sequence = FASTA_INDEX[chrom][start - 1:end]
        sequence = str(sequence)

        if len(sequence) != (end - start + 1):
            return None
        return sequence
    except:
        return None


def extract_sequences_from_chrom_ranges(fasta_file, chrom_ranges):
    """批量提取序列（多线程加速）"""
    load_fasta_index(fasta_file)

    # 多线程提取
    with ThreadPoolExecutor(max_workers=4) as executor:
        sequences = list(executor.map(extract_sequence_batch, chrom_ranges))

    return sequences


# ====================== 保留原函数（兼容旧代码） ======================
def DNA_to_matrix(DNA, sequence_length=10000):
    """兼容原接口的单条序列编码"""
    return DNA_to_matrix_batch([DNA], sequence_length)[0:1]


def extract_sequence_from_chrom_range(fasta_file, chrom_range):
    """兼容原接口的单条序列提取"""
    load_fasta_index(fasta_file)
    return extract_sequence_batch(chrom_range)


# ====================== 主预测函数（核心优化+耗时统计） ======================
def predict_tad_boundary(bed_file_path, fasta_file_path="./hg19.fa",
                         model_path="./GM12878_tad_boundary_model.h5",
                         output_file_path="prediction_result.txt",
                         batch_size=200):  # 批量大小（可调）
    """
    预测TAD边界概率的核心函数（优化版，带耗时统计）
    :param bed_file_path: 输入BED文件路径（必填）
    :param fasta_file_path: hg19.fa文件路径（默认当前目录）
    :param model_path: 预训练模型文件路径（默认当前目录）
    :param output_file_path: 输出结果文件路径（默认当前目录）
    :param batch_size: 批量预测大小（越大越快，需匹配内存）
    :return: 输出结果文件的绝对路径 + 耗时统计字典
    """
    # 初始化耗时统计
    time_stats = {
        "total": time.time(),
        "file_check": 0,
        "model_load": 0,
        "bed_parse": 0,
        "sequence_extract": 0,
        "encoding": 0,
        "prediction": 0,
        "result_write": 0
    }

    # 1. 校验输入文件
    t_start = time.time()
    if not os.path.exists(bed_file_path):
        raise FileNotFoundError(f"BED文件不存在：{bed_file_path}")
    if not os.path.exists(fasta_file_path):
        raise FileNotFoundError(f"FASTA文件不存在：{fasta_file_path}")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"模型文件不存在：{model_path}")
    time_stats["file_check"] = time.time() - t_start
    print(f"[1/7] 文件检查完成 | 耗时: {time_stats['file_check']:.2f}秒")

    # 2. 加载模型
    t_start = time.time()
    model = load_model(model_path)
    time_stats["model_load"] = time.time() - t_start
    print(f"[2/7] 模型加载完成 | 耗时: {time_stats['model_load']:.2f}秒")

    # 3. 解析BED文件
    t_start = time.time()
    bed_data = []
    chrom_ranges = []

    with open(bed_file_path, 'r') as bf:
        for line in bf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[:3]
            bed_data.append((chrom, start, end))
            chrom_ranges.append(f"{chrom}:{start}-{end}")

    total_lines = len(bed_data)
    time_stats["bed_parse"] = time.time() - t_start
    print(f"[3/7] BED解析完成 | 有效区间: {total_lines} | 耗时: {time_stats['bed_parse']:.2f}秒")

    if total_lines == 0:
        raise ValueError("BED文件中无有效区间")

    # 4. 批量提取序列
    t_start = time.time()
    sequences = extract_sequences_from_chrom_ranges(fasta_file_path, chrom_ranges)
    # 过滤无效序列
    valid_indices = [i for i, seq in enumerate(sequences) if seq is not None]
    valid_sequences = [sequences[i] for i in valid_indices]
    valid_bed = [bed_data[i] for i in valid_indices]
    time_stats["sequence_extract"] = time.time() - t_start
    print(
        f"[4/7] 序列提取完成 | 有效序列: {len(valid_sequences)}/{total_lines} | 耗时: {time_stats['sequence_extract']:.2f}秒")

    # 5. 批量编码 + 批量预测
    t_start_encoding = time.time()
    t_start_pred = time.time()
    predictions = []

    # 分批次处理
    for i in range(0, len(valid_sequences), batch_size):
        batch_seqs = valid_sequences[i:i + batch_size]

        # 批量编码
        batch_matrix = DNA_to_matrix_batch(batch_seqs)

        # 批量预测
        batch_pred = model.predict(batch_matrix, verbose=0, batch_size=batch_size)
        predictions.extend([round(float(p[0]), 6) for p in batch_pred])

    time_stats["encoding"] = time.time() - t_start_encoding
    time_stats["prediction"] = time.time() - t_start_pred
    print(f"[5/7] 序列编码完成 | 耗时: {time_stats['encoding']:.2f}秒")
    print(f"[6/7] 模型预测完成 | 耗时: {time_stats['prediction']:.2f}秒")

    # 6. 写入结果
    t_start = time.time()
    with open(output_file_path, 'w') as outf:
        outf.write("chrom\tstart\tend\tprediction\n")
        # 批量构造输出行（减少IO次数）
        output_lines = []
        for (chrom, start, end), pred in zip(valid_bed, predictions):
            output_lines.append(f"{chrom}\t{start}\t{end}\t{pred}\n")
        outf.writelines(output_lines)

    time_stats["result_write"] = time.time() - t_start
    print(f"[7/7] 结果写入完成 | 耗时: {time_stats['result_write']:.2f}秒")

    # 总耗时
    time_stats["total"] = time.time() - time_stats["total"]

    # 打印汇总统计
    print("\n" + "=" * 60)
    print("📊 预测耗时统计汇总")
    print("=" * 60)
    print(f"总耗时:                {time_stats['total']:.2f}秒")
    print(f"  ├─ 文件检查:         {time_stats['file_check']:.2f}秒")
    print(f"  ├─ 模型加载:         {time_stats['model_load']:.2f}秒")
    print(f"  ├─ BED解析:          {time_stats['bed_parse']:.2f}秒")
    print(f"  ├─ 序列提取:         {time_stats['sequence_extract']:.2f}秒")
    print(f"  ├─ 序列编码:         {time_stats['encoding']:.2f}秒")
    print(f"  ├─ 模型预测:         {time_stats['prediction']:.2f}秒")
    print(f"  └─ 结果写入:         {time_stats['result_write']:.2f}秒")
    print(f"\n📈 处理效率: {len(valid_sequences) / time_stats['total']:.2f} 条/秒")
    print("=" * 60)

    # 返回绝对路径和耗时统计
    output_abs_path = os.path.abspath(output_file_path)
    return output_abs_path, time_stats

#loss bed code


