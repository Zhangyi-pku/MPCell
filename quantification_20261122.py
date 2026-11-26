#!/usr/bin/env python
# coding: utf-8

import sys
import json
import os
import numpy as np
import pandas as pd
from scipy.special import comb

# 强制实时输出
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
else:
    import functools
    print = functools.partial(print, flush=True)

# ------------------------------------------------------------
# 辅助函数：ppm 转绝对质量容差
# ------------------------------------------------------------
def ppm2dalton(mz, ppm, res_define_at=200):
    """计算给定 m/z 的 ppm 绝对容差（与 MetCell 逻辑一致）"""
    if mz >= res_define_at:
        return mz * ppm
    else:
        return res_define_at * ppm

# ------------------------------------------------------------
# 模拟理论同位素强度比（使用烷烃近似）
# ------------------------------------------------------------
def _simulate_theo_isotope(mz, isotope_max_num=4):
    """
    根据 m/z 近似估计碳原子数（假设为烷烃 CnH2n+2），
    返回各同位素等级的相对强度（单同位素为 1）。
    """
    n_carbon = max(1, int((mz - 2*1.0078) / 14.0156))
    R = 0.011  # 13C 自然丰度，忽略氘
    intensities = [1.0]
    for i in range(1, isotope_max_num):
        intensity = comb(n_carbon, i) * (R ** i)
        intensities.append(intensity)
    return np.array(intensities)

# ------------------------------------------------------------
# 主函数：同位素指认
# ------------------------------------------------------------
def annotate_isotope(df,
                     mz_tol_ppm=10e-6,
                     isotope_delta=1.003355,
                     isotope_max_num=4,
                     isotope_int_ratio_check=True,
                     isotope_int_ratio_cutoff=500.0,
                     res_define_at=200):
    df = df.copy()
    df['base_peak'] = pd.NA
    df['isotope'] = pd.NA
    df = df.sort_values('mean_intensity', ascending=False).reset_index(drop=True)
    df.insert(0, 'name', df['mz'].tolist())

    unmarked = df['base_peak'].isna()
    total = unmarked.sum()
    processed = 0
    print(f"Starting isotope annotation, {total} candidate peak in total", flush=True)

    for i, row in df[unmarked].iterrows():
        processed += 1
        if processed % max(1, total // 20) == 0 or processed == total:
            percent = processed / total * 100
            print(f"Schedule of isotope annotation: {processed}/{total} ({percent:.1f}%)", flush=True)

        if not pd.isna(row['base_peak']):
            continue
        base_name = row['name']
        base_mz = row['mz']
        base_intensity = row['mean_intensity']

        # 若 base 本身强度为 0（未匹配峰），直接标 [M] 并跳过候选搜索，
        # 避免 0/0 或 0 强度候选干扰同位素判定
        if base_intensity <= 0:
            df.at[i, 'base_peak'] = base_name
            df.at[i, 'isotope'] = '[M]'
            continue

        mz_isotope_list = [base_mz + i * isotope_delta for i in range(isotope_max_num)]
        mz_windows = []
        for iso_mz in mz_isotope_list:
            tol = ppm2dalton(iso_mz, mz_tol_ppm, res_define_at)
            mz_windows.append((iso_mz - tol, iso_mz + tol))

        candidates = []
        for j, cand in df[unmarked].iterrows():
            if cand['name'] == base_name:
                continue
            # [FIX] 跳过 0 强度候选（未匹配峰），避免其被误标为 [M+k]
            if cand['mean_intensity'] <= 0:
                continue
            for k, (low, high) in enumerate(mz_windows):
                if low <= cand['mz'] <= high:
                    candidates.append({
                        'index': j,
                        'name': cand['name'],
                        'mz': cand['mz'],
                        'mean_intensity': cand['mean_intensity'],
                        'iso_level': k
                    })
                    break

        if not candidates:
            df.at[i, 'base_peak'] = base_name
            df.at[i, 'isotope'] = '[M]'
            continue

        theo_intensities = _simulate_theo_isotope(base_mz, isotope_max_num)
        accepted = []
        for cand in candidates:
            iso_lvl = cand['iso_level']
            obs_ratio = cand['mean_intensity'] / base_intensity
            theo_ratio = theo_intensities[iso_lvl]
            if isotope_int_ratio_check:
                if theo_ratio == 0:
                    continue
                error = abs(obs_ratio - theo_ratio) / theo_ratio * 100
                if error <= isotope_int_ratio_cutoff:
                    accepted.append(cand)
            else:
                accepted.append(cand)

        best_by_level = {}
        for cand in accepted:
            lvl = cand['iso_level']
            if lvl not in best_by_level:
                best_by_level[lvl] = cand
            else:
                curr_err = abs(best_by_level[lvl]['mz'] - mz_isotope_list[lvl])
                new_err = abs(cand['mz'] - mz_isotope_list[lvl])
                if new_err < curr_err:
                    best_by_level[lvl] = cand

        df.at[i, 'base_peak'] = base_name
        df.at[i, 'isotope'] = '[M]'
        for lvl, cand in best_by_level.items():
            idx_cand = cand['index']
            df.at[idx_cand, 'base_peak'] = base_name
            df.at[idx_cand, 'isotope'] = f'[M+{lvl}]'

    print("Finish", flush=True)
    return df

def main():
    if len(sys.argv) < 2:
        print("Error! No paras_file path specified", file=sys.stderr, flush=True)
        sys.exit(1)

    param_file = sys.argv[1]
    try:
        with open(param_file, 'r', encoding='utf-8') as f:
            params = json.load(f)
    except Exception as e:
        print(f"Failed to read paras_file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # 获取基本参数
    mz_list_path = params.get('mz_list', '')
    output_dir = params.get('output_dir', '')
    ppm = params.get('ppm', 10e-6)
    print(f"(Read from paras_file) ppm = {params.get('ppm')}", flush=True)

    if not mz_list_path:
        print("Error! No peak_list path specified", file=sys.stderr, flush=True)
        sys.exit(1)
    if not output_dir:
        print("Error! No output path specified", file=sys.stderr, flush=True)
        sys.exit(1)

    # 构建路径
    path_p = os.path.join(output_dir, 'grouping', 'RT_del')
    c_peak_file = os.path.join(output_dir, 'c_peak.npy')
    blank_aver_file = os.path.join(path_p, 'c_peak_blank_aver_intensity.npy')

    # 检查必要文件是否存在
    missing = []
    if not os.path.exists(c_peak_file):
        missing.append(c_peak_file)
    if not os.path.exists(blank_aver_file):
        missing.append(blank_aver_file)
    if missing:
        print("Error! Missed files listed: ", file=sys.stderr, flush=True)
        for f in missing:
            print(f"  {f}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ------------------------------------------------------------
    # 1. 读取用户提供的参照峰列表（Excel 文件，仅 m/z）
    # ------------------------------------------------------------
    try:
        df_user = pd.read_excel(mz_list_path)
        if 'm/z' not in df_user.columns and 'mz' not in df_user.columns:
            print("Error! Excel must include 'm/z' or 'mz' columns", file=sys.stderr, flush=True)
            sys.exit(1)
        mz_col = 'm/z' if 'm/z' in df_user.columns else 'mz'
        ref_peaks = df_user[mz_col].tolist()
        print(f"Read {len(ref_peaks)} peaks (m/z)", flush=True)
    except Exception as e:
        print(f"Failed to read provided file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ------------------------------------------------------------
    # 2. 读取 c_peak.npy（no_tims 结构：[mz, rt_list, int_list]）
    # ------------------------------------------------------------
    try:
        c_peak = np.load(c_peak_file, allow_pickle=True).tolist()
        print(f"c_peak length {len(c_peak)}", flush=True)
    except Exception as e:
        print(f"Failed to read c_peak.npy: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ------------------------------------------------------------
    # 3. 读取背景平均强度
    # ------------------------------------------------------------
    try:
        blank_aver_intensity = np.load(blank_aver_file)
        print(f"Blank_aver_peak length: {len(blank_aver_intensity)}", flush=True)
    except Exception as e:
        print(f"Failed to read c_peak_blank_aver_intensity.npy: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ------------------------------------------------------------
    # 4. 构建强度矩阵
    # ------------------------------------------------------------
    intensity_matrix = []
    matched_count = 0
    unmatched = []

    mzppm = ppm

    # [FIX] 预先求每行应有的长度（= pulse 数），用于未匹配行占位，
    # 避免 pd.DataFrame 对 ragged list（含空列表）报 ValueError。
    n_pulses_total = len(c_peak[0][2]) if len(c_peak) > 0 else 0
    print(f"Each intensity row length = {n_pulses_total} (number of pulses)", flush=True)

    for idx_ref, mz_ref in enumerate(ref_peaks):
        mz_tol = mz_ref * mzppm
        found = False
        for idx_peak, peak_info in enumerate(c_peak):
            # no_tims 结构：[mz, rt_list, int_list]
            mz_peak = peak_info[0]

            if abs(mz_ref - mz_peak) <= mz_tol:
                intensity_list = peak_info[2]
                # 处理 None（worker 中当 pulse 内无有效 RT 时可能为 None）
                intensity_list = [0 if x is None else x for x in intensity_list]
                blank_val = blank_aver_intensity[idx_peak]
                peak_intensity = np.array(intensity_list, dtype=float) - blank_val
                peak_intensity[peak_intensity < 0] = 0
                intensity_matrix.append(peak_intensity.tolist())
                found = True
                matched_count += 1
                break
        if not found:
            print(f"Error! No matched peak m/z = {mz_ref}", flush=True)
            unmatched.append(mz_ref)
            # [FIX] 用等长全 0 占位，保证 intensity_matrix 规整
            intensity_matrix.append([0.0] * n_pulses_total)

    print(f" {matched_count} / {len(ref_peaks)} peak matched successfully", flush=True)

    # ------------------------------------------------------------
    # 5. 保存强度矩阵 Excel
    # ------------------------------------------------------------
    df_intensity = pd.DataFrame(intensity_matrix)
    mz_values = [float(mz) for mz in ref_peaks]
    df_intensity.insert(0, 'mz', mz_values)
    output_excel = os.path.join(output_dir, 'intensity_matrix.xlsx')
    try:
        df_intensity.to_excel(output_excel)
        print(f"Intensity matrix saved in: {output_excel}", flush=True)
    except Exception as e:
        print(f"Failed to save excel file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # 6. 如果有未匹配的峰，写入单独文件
    if unmatched:
        unmatched_path = os.path.join(output_dir, 'unmatched_mz.txt')
        with open(unmatched_path, 'w') as f:
            for mz in unmatched:
                f.write(f"{mz}\n")
        print(f"Unmatched peak saved in: {unmatched_path}", flush=True)

    # ------------------------------------------------------------
    # 7. 同位素指认
    # ------------------------------------------------------------
    print("Starting isotope annotation", flush=True)
    ref_mean_intensity = []
    for row in intensity_matrix:
        if len(row) == 0:
            ref_mean_intensity.append(0.0)
            continue
        arr = np.array(row, dtype=float)
        nonzero = arr[arr != 0]
        if len(nonzero) != 0:
            ref_mean_intensity.append(float(np.mean(nonzero)))
        else:
            ref_mean_intensity.append(0.0)
    df_intensity['mean_intensity'] = ref_mean_intensity

    iso_result = annotate_isotope(df_intensity, mz_tol_ppm=ppm)
    iso_output_excel = os.path.join(output_dir, 'intensity_matrix_iso_anno.xlsx')
    try:
        iso_result.to_excel(iso_output_excel)
        print(f"Intensity matrix with isotope annotation saved in: {iso_output_excel}", flush=True)
    except Exception as e:
        print(f"Failed to save excel file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ---------- 筛选只保留 [M] 同位素峰（单同位素峰） ----------
    if 'isotope' in iso_result.columns:
        df_mono = iso_result[iso_result['isotope'] == '[M]'].copy()
        df_mono_final = df_mono.drop(columns=['name','mean_intensity','base_peak','isotope'])
        mono_output_excel = os.path.join(output_dir, 'intensity_matrix_iso_anno_M0.xlsx')
        try:
            df_mono_final.to_excel(mono_output_excel, index=False)
            print(f"Intensity matrix with only [M] isotopes saved in: {mono_output_excel}", flush=True)
        except Exception as e:
            print(f"Failed to save mono isotope file: {e}", file=sys.stderr, flush=True)
            sys.exit(1)
    else:
        print("Warning: 'isotope' column not found in iso_result, skipping mono isotope extraction.", flush=True)

    print("Quantification finished", flush=True)

if __name__ == "__main__":
    main()