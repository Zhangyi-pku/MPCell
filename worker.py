#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sys
import json
import os
import bisect
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # 防止在没有显示设备时报错
import matplotlib.pyplot as plt
from pyteomics import mzml
from scipy.signal import find_peaks
from statsmodels.nonparametric.smoothers_lowess import lowess
import heapq
import copy
import time
import decimal
import bisect


# In[ ]:


# 强制 stdout/stderr 行缓冲，确保实时输出
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
else:
    import functools
    print = functools.partial(print, flush=True)

# ------------------------------------------------------------
# 模块 1: 读取 mzML 文件，提取一级质谱离子列表
# ------------------------------------------------------------
def read_mzml_and_extract_ions(input_file, peak_intensity_thres, rt_range_max):
    print(f"[1/14] Read mzML file: {input_file}", flush=True)
    reader = mzml.read(input_file)
    rt = []
    rt_to_ions = {}
    for spectrum in reader:
        if spectrum['ms level'] == 1:
            rt0 = float(spectrum['scanList']['scan'][0]['scan start time'])
            if len(rt) == 0 or rt0 != rt[-1]:
                rt.append(rt0)
            mz_arr = np.array(spectrum['m/z array'])
            int_arr = np.array(spectrum['intensity array'])
            mask = int_arr > peak_intensity_thres
            if np.any(mask):
                mz_arr = mz_arr[mask]
                int_arr = int_arr[mask]
                order = np.argsort(mz_arr)
                rt_to_ions[rt0] = (mz_arr[order], int_arr[order])
            if rt0 > rt_range_max:
                break
    print(f"     Read {len(rt)} scans (MS1 spec)，effective spec: {len(rt_to_ions)}", flush=True)
    print("[1/14] Finish", flush=True)
    return rt, rt_to_ions

# ------------------------------------------------------------
# 模块 2: 根据 rt_range 筛选保留时间范围
# ------------------------------------------------------------
def filter_rt_range(rt, rt_range_min, rt_range_max):
    print(f"[2/14] RT range [{rt_range_min}, {rt_range_max}]", flush=True)
    rt = np.array(rt)
    mask = (rt > rt_range_min) & (rt < rt_range_max)
    filtered_rt = rt[mask]
    filtered_index = np.arange(len(filtered_rt))
    print(f"     {len(filtered_rt)} time points reserved", flush=True)
    print("[2/14] Finish", flush=True)
    return filtered_rt, filtered_index

# ------------------------------------------------------------
# 模块 3: 提取指定 marker 的强度曲线
# ------------------------------------------------------------
def extract_marker_intensity(rt, rt_to_ions, marker, ppm):
    print(f"[3/14] Extract marker m/z = {marker} (ppm={ppm}) EIC", flush=True)
    marker_intensity = np.zeros(len(rt))
    detected_count = 0
    for idx, rt_val in enumerate(rt):
        if rt_val not in rt_to_ions:
            continuecell_marker
        mz_arr, int_arr = rt_to_ions[rt_val]
        tol = marker * ppm
        mask = abs(mz_arr - marker) < tol
        if np.any(mask):
            marker_intensity[idx] = int_arr[mask][0]
            detected_count += 1
    print("[3/14] Finish", flush=True)
    return marker_intensity

# ------------------------------------------------------------
# 模块 4: 基于 marker 强度检测细胞脉冲
# ------------------------------------------------------------
def detect_cell_pulses(rt, index, marker_intensity, intensity_thres, cell_spec):
    print(f"[4/14] Cell pulse detection (intensity threshold={intensity_thres}, cell_spec={cell_spec})", flush=True)
    mask = marker_intensity > intensity_thres
    mask_idx = np.where(mask)[0]
    if len(mask_idx) == 0:
        print("     ⚠️ No cell pulse detected, processing terminated", flush=True)
        print("[4/14] Finish (no cell pulses)", flush=True)
        return [], [], [], [], []
    cell_rt = rt[mask_idx]
    cell_index = index[mask_idx]
    cell_intensity = marker_intensity[mask_idx]

    diff_idx = np.diff(cell_index)
    break_points = np.where(diff_idx != 1)[0] + 1
    groups = np.split(mask_idx, break_points)

    if cell_spec % 2 == 0:
        spec_range = range(int(-cell_spec/2), int(cell_spec/2))
    else:
        spec_range = range(int(-cell_spec/2), int(cell_spec/2) + 1)

    cell_pulse = []
    cell_marker_vals = []
    for each_group in groups:
        group_intensity = marker_intensity[each_group]
        local_max_idx = np.argmax(group_intensity)
        global_max_idx = each_group[local_max_idx]
        start = global_max_idx + spec_range[0]
        end = global_max_idx + spec_range[-1]
        if start > 0 and end < len(rt):
            pulse_rts = rt[global_max_idx + np.array(spec_range)].tolist()
            cell_pulse.append(pulse_rts)
            pulse_marker = marker_intensity[global_max_idx + np.array(spec_range)].tolist()
            cell_marker_vals.append(pulse_marker)
    print(f"     {len(cell_pulse)} cell pulses detected", flush=True)
    print("[4/14] Finish", flush=True)
    return cell_pulse, cell_marker_vals, cell_rt, cell_index, cell_intensity

# ------------------------------------------------------------
# 模块 5: 迭代聚类，筛选特征峰 m/z 列表
# ------------------------------------------------------------
def cluster_features(rt_list, rt_to_ions, ppm, cell_pulse):
    print(f"[5/14] peak grouping (ppm={ppm}), {len(cell_pulse)} cell pulses used", flush=True)
    mz_list = []
    for rt_val in rt_list:
        mz_arr, int_arr = rt_to_ions[rt_val]
        for idx, mz in enumerate(mz_arr):
            mz_list.append((mz, int_arr[idx]))
    mz_list.sort(key=lambda x: x[1], reverse=True)
    mz_list = np.array(mz_list)
    print(f"     Original mz_list length: {len(mz_list)}", flush=True)

    cycle = 0
    while True:
        clusters = []
        ref_mzs = []
        for (mz, intensity) in mz_list:
            pos = bisect.bisect_left(ref_mzs, mz)
            found = False
            candidate_int = -1
            candidate_idx = -1
            for idx in (pos, pos-1):
                if 0 <= idx < len(clusters):
                    ref = clusters[idx]['ref_mz']
                    if abs(mz - ref) <= ref * ppm:
                        max_int = clusters[idx]['points'][0][1]
                        if max_int > candidate_int:
                            candidate_idx = idx
                            candidate_int = max_int
                        found = True
            if found:
                clusters[candidate_idx]['points'].append((mz, intensity))
            else:
                new_cluster = {'ref_mz': mz, 'points': [(mz, intensity)]}
                clusters.insert(pos, new_cluster)
                ref_mzs.insert(pos, mz)

        percent = 0.1 if cycle == 0 else 0
        mz_list_final = []
        for cluster in clusters:
            points = np.array(cluster['points'])
            if len(points) > percent * len(cell_pulse):
                plus_sum = np.sum(points[:, 0] * points[:, 1])
                intensity_sum = np.sum(points[:, 1])
                mz_tar = plus_sum / intensity_sum
                mz_list_final.append((mz_tar, intensity_sum))
        mz_list_final.sort(key=lambda x: x[1], reverse=True)
        cycle += 1
        print(f"     iteration {cycle}: group={len(clusters)}, mz_list length={len(mz_list_final)}", flush=True)
        if len(mz_list) == len(mz_list_final):
            break
        mz_list = np.array(mz_list_final)

    mz_tar = np.array([x[0] for x in mz_list_final])
    print(f"     Final mz_list length: {len(mz_tar)} ", flush=True)
    print("[5/14] Finish", flush=True)
    return mz_tar

# ------------------------------------------------------------
# 模块 6: 构建 c_peak 的 EIC 强度矩阵
# ------------------------------------------------------------
def extract_eic_matrix(rt, rt_to_ions, c_peak_mz, ppm):
    print(f"[6/14] Constrction of EIC intensity matrix (peak length={len(c_peak_mz)}, time point={len(rt)})", flush=True)
    n_times = len(rt)
    n_peaks = len(c_peak_mz)
    sorted_idx = np.argsort(c_peak_mz)
    c_peak_mz_sorted = c_peak_mz[sorted_idx]
    orig_order = np.argsort(sorted_idx)

    time_intensity = np.zeros((n_peaks, n_times), dtype=np.float32)
    for t_idx, rt0 in enumerate(rt):
        if rt0 not in rt_to_ions:
            continue
        mz_arr, int_arr = rt_to_ions[rt0]
        for peak_idx in range(n_peaks):
            mz_target = c_peak_mz_sorted[peak_idx]
            tol = mz_target * ppm
            left = bisect.bisect_left(mz_arr, mz_target - tol)
            right = bisect.bisect_right(mz_arr, mz_target + tol)
            if left < right:
                time_intensity[peak_idx, t_idx] = int_arr[left:right][0]
    time_intensity_orig = time_intensity[orig_order, :]
    print(f"     EIC shape: {time_intensity_orig.shape}", flush=True)
    print("[6/14] Finish", flush=True)
    return time_intensity_orig

# ------------------------------------------------------------
# 模块 7: 对每个细胞脉冲，提取该峰的最大强度及对应 RT
# ------------------------------------------------------------
def extract_cell_peak_max(rt, cell_pulse, eic_matrix, c_peak_mz):
    print(f"[7/14] Max intensity extraction in cell pulses({len(cell_pulse)} cell pulses, {len(c_peak_mz)} peaks)", flush=True)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(rt)}
    n_peaks = len(c_peak_mz)
    c_peak_max = []
    for each in range(n_peaks):
        intensities = eic_matrix[each, :]
        rt_max_list = []
        intensity_max_list = []
        for pulse in cell_pulse:
            max_intensity = -1
            max_rt = None
            for rt0 in pulse:
                idx = rt_to_index.get(rt0)
                if idx is not None:
                    val = intensities[idx]
                    if val > max_intensity:
                        max_intensity = val
                        max_rt = rt0
            if max_rt is not None:
                rt_max_list.append(max_rt)
                intensity_max_list.append(max_intensity)
            else:
                rt_max_list.append(None)
                intensity_max_list.append(None)
        c_peak_max.append([c_peak_mz[each], rt_max_list, intensity_max_list])
    print("[7/14] Finish", flush=True)
    return np.array(c_peak_max, dtype=object)

# ------------------------------------------------------------
# 模块 8: 将时间轴分为三组，并提取子集数据
# ------------------------------------------------------------
def create_groups_and_subset(output_dir, rt, cell_pulse, c_peak, eic_matrix,
                             bg_marker, bg_intensity_max, ppm, cell_spec_lim, marker):
    print(f"[8/14] RT subgrouping and subsets extraction (RT_del)", flush=True)
    # 确保 cell_pulse 和 c_peak 为 numpy 数组（以便布尔索引）
    cell_pulse = np.array(cell_pulse, dtype=object)
    c_peak = np.array(c_peak, dtype=object)

    rt_list = rt.tolist()
    rt_del = int(rt_list[-1] // 4)
    if rt_del > 5:
        rt_del = 5
    RT_5min = [r for r in rt_list if r > rt_del]
    RT_interval = int(len(RT_5min) / 3)
    RT_part1 = RT_5min[0:RT_interval]
    RT_part2 = RT_5min[RT_interval:RT_interval*2]
    RT_part3 = RT_5min[RT_interval*2:len(RT_5min)]

    path_group = os.path.join(output_dir, 'grouping')
    os.makedirs(path_group, exist_ok=True)
    np.save(os.path.join(path_group, 'RT_part1.npy'), RT_part1)
    np.save(os.path.join(path_group, 'RT_part2.npy'), RT_part2)
    np.save(os.path.join(path_group, 'RT_part3.npy'), RT_part3)
    np.save(os.path.join(path_group, 'RT_del.npy'), RT_5min)

    part = 'RT_del'
    path_p = os.path.join(path_group, part)
    os.makedirs(path_p, exist_ok=True)
    RT_sub = np.load(os.path.join(path_group, part+'.npy')).tolist()
    np.save(os.path.join(path_p, 'RT.npy'), RT_sub)

    rt_min, rt_max = RT_sub[0], RT_sub[-1]
    start_times = np.array([p[0] for p in cell_pulse])
    end_times = np.array([p[-1] for p in cell_pulse])
    mask = (start_times >= rt_min) & (end_times <= rt_max)

    cell_pulse_indices = np.where(mask)[0]
    cell_pulse_p = cell_pulse[mask]

    if len(cell_pulse_indices) == 0:
        raise ValueError("该分组中没有细胞脉冲，请检查 RT 范围")

    start_cell_ind = cell_pulse_indices[0]
    end_cell_ind = cell_pulse_indices[-1]

    c_peak_p = []
    for i in range(len(c_peak)):
        c_peak_p.append([c_peak[i][0]] + [arr[start_cell_ind:end_cell_ind+1] for arr in c_peak[i][1:]])

    rt_mask = (rt >= rt_min) & (rt <= rt_max)
    RT_indices = np.where(rt_mask)[0]
    start_rt_ind = RT_indices[0]
    end_rt_ind = RT_indices[-1]

    eic_matrix_p = eic_matrix[:, start_rt_ind:end_rt_ind+1]

    np.save(os.path.join(path_p, 'cell_pulse.npy'), np.array(cell_pulse_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak.npy'), np.array(c_peak_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak_eic_intensity.npy'), eic_matrix_p)

    print(f"     subgroup: RT range [{rt_min:.2f}, {rt_max:.2f}], cell pulse={len(cell_pulse_p)}, peak={len(c_peak_p)}", flush=True)
    print("[8/14] Finish", flush=True)
    return path_p, cell_pulse_p, c_peak_p, eic_matrix_p, RT_sub

# ------------------------------------------------------------
# 模块 9: 计算细胞脉冲内峰值平均强度和 RSD
# ------------------------------------------------------------
def compute_average_intensity_and_rsd(c_peak_p):
    print(f"[9/14] Peak information in cell pulses", flush=True)
    n_peaks = len(c_peak_p)
    aver_intensity = np.zeros(n_peaks)
    rsd = np.zeros(n_peaks)
    for i in range(n_peaks):
        intensity = c_peak_p[i][2]
        if not isinstance(intensity, np.ndarray):
            intensity = np.array(intensity)
        non_zero = intensity[intensity != 0]
        if len(non_zero) > 0:
            sorted_vals = np.sort(non_zero)
            trim_len = max(1, int(0.1 * len(sorted_vals)))
            trimmed = sorted_vals[trim_len:-trim_len] if len(sorted_vals) > 2*trim_len else sorted_vals
            if len(trimmed) > 0:
                mean_val = np.mean(trimmed)
                aver_intensity[i] = mean_val
                if mean_val != 0:
                    rsd[i] = np.std(trimmed, ddof=1) / mean_val
    print(f"     Finish, averaged intensity range: [{np.min(aver_intensity):.2e}, {np.max(aver_intensity):.2e}], averaged RSD: {np.mean(rsd):.4f}", flush=True)
    print("[9/14] Finish", flush=True)
    return aver_intensity, rsd

# ------------------------------------------------------------
# 模块 10: 计算背景（非细胞区域）的平均强度、标准差、RSD
# ------------------------------------------------------------
def compute_background_stats(path_p, bg_marker, bg_intensity_max, ppm, c_peak_mz_all, eic_matrix_p, cell_pulse_p, RT_sub):
    print(f"[10/14] Peak information in background (bg_marker={bg_marker}, bg_intensity_max={bg_intensity_max})", flush=True)
    # 找出所有背景时间点
    cell_pulse_rt = set()
    for pulse in cell_pulse_p:
        cell_pulse_rt.update(pulse)
    cell_pulse_rt = sorted(cell_pulse_rt)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(RT_sub)}
    cell_pulse_idx = [rt_to_index[rt0] for rt0 in cell_pulse_rt if rt0 in rt_to_index]
    bg_nocell_idx = [idx for idx in range(len(RT_sub)) if idx not in cell_pulse_idx]
    print(f"     no-cell BG time point: {len(bg_nocell_idx)} / {len(RT_sub)}", flush=True)

    c_peak_blank_index = None
    for count, (bg_mz, bg_max) in enumerate(zip(bg_marker, bg_intensity_max)):
        mask = np.abs(c_peak_mz_all - bg_mz) / bg_mz < ppm
        peak_indices = np.where(mask)[0]
        if len(peak_indices) == 0:
            raise ValueError(f"No peak matches background marker {bg_mz}")
        eic_intensity_marker = eic_matrix_p[peak_indices[0]]
        bg_intensity_idx = np.where(eic_intensity_marker <= bg_max)[0]
        if count == 0:
            c_peak_blank_index = set(bg_intensity_idx)
        else:
            c_peak_blank_index.intersection_update(bg_intensity_idx)
        print(f"     BG marker {bg_mz} filtered BG time point: {len(c_peak_blank_index)}", flush=True)

    c_peak_blank_index = np.array(list(c_peak_blank_index))
    c_peak_blank_index = np.intersect1d(bg_nocell_idx, c_peak_blank_index)
    print(f"     Final BG time point: {len(c_peak_blank_index)}", flush=True)

    n_peaks = eic_matrix_p.shape[0]
    blank_aver = np.zeros(n_peaks)
    blank_std = np.zeros(n_peaks)
    blank_rsd = np.zeros(n_peaks)
    for i in range(n_peaks):
        blank_intensity = eic_matrix_p[i, c_peak_blank_index]
        non_zero = blank_intensity[blank_intensity != 0]
        if len(non_zero) > min(50,max(10, 0.01*len(c_peak_blank_index))):
            sorted_vals = np.sort(non_zero)
            trim_len = int(0.1*len(sorted_vals))
            trimmed = sorted_vals[trim_len:-trim_len] if len(sorted_vals) > 2*trim_len else sorted_vals
        else:
            sorted_vals = np.sort(blank_intensity)
            trim_len = int(0.1*len(sorted_vals))
            trimmed = sorted_vals[trim_len:-trim_len] if len(sorted_vals) > 2*trim_len else sorted_vals
        mean_val = np.mean(trimmed)
        blank_aver[i] = mean_val
        blank_std[i] = np.std(trimmed, ddof=1)
        if mean_val != 0:
            blank_rsd[i] = blank_std[i] / mean_val
    print(f"     avergaed intensity range: [{np.min(blank_aver):.2e}, {np.max(blank_aver):.2e}]", flush=True)
    print("[10/14] Finish", flush=True)
    return blank_aver, blank_std, blank_rsd, c_peak_blank_index

# ------------------------------------------------------------
# 模块 11: 计算 cell_occur（细胞中出现次数）
# ------------------------------------------------------------
def compute_cell_occurrence(eic_matrix_p, cell_pulse_p, blank_aver, blank_std, RT_sub, cell_spec_lim):
    print(f"[11/14] Cell occurrence (cell_spec_lim={cell_spec_lim})", flush=True)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(RT_sub)}
    n_peaks = eic_matrix_p.shape[0]
    cell_occur = np.zeros(n_peaks, dtype=np.int64)
    for each in range(n_peaks):
        intensities = eic_matrix_p[each, :]
        blank_val = blank_aver[each]
        blank_sd = blank_std[each]
        count = 0
        for pulse in cell_pulse_p:
            intensity_list = []
            for rt0 in pulse:
                idx = rt_to_index.get(rt0)
                if idx is not None:
                    intensity_list.append(intensities[idx])
            idx1 = np.where(np.array(intensity_list) > blank_val + 3 * blank_sd)[0]
            if len(idx1) > 0:
                diff_idx = np.diff(idx1)
                if len(np.where(diff_idx == 1)[0]) >= cell_spec_lim - 1:
                    count += 1
        cell_occur[each] = count
    print(f"     Cell occurrence summary: min={np.min(cell_occur)}, max={np.max(cell_occur)}, mean={np.mean(cell_occur):.2f}", flush=True)
    print("[11/14] Finish", flush=True)
    return cell_occur

# ------------------------------------------------------------
# 模块 12: 计算 c_peak_pulse（脉冲数）
# ------------------------------------------------------------
def compute_peak_pulse_count(eic_matrix_p, blank_aver, blank_std):
    print(f"[12/14] Peak pulse calculation", flush=True)
    n_peaks = eic_matrix_p.shape[0]
    pulse1 = np.zeros(n_peaks, dtype=np.int64)
    pulse2 = np.zeros(n_peaks, dtype=np.int64)
    for i in range(n_peaks):
        intensity = eic_matrix_p[i]
        thresh = blank_aver[i] + 3 * blank_std[i]
        pos = np.where(intensity > thresh)[0]
        if len(pos) > 0:
            diff = np.diff(pos)
            split_indices = np.where(diff != 1)[0] + 1
            groups = np.split(pos, split_indices)
            pulse1[i] = len(groups)
        pos0 = np.where(intensity > 0)[0]
        if len(pos0) > 0:
            diff0 = np.diff(pos0)
            split_indices0 = np.where(diff0 != 1)[0] + 1
            groups0 = np.split(pos0, split_indices0)
            pulse2[i] = len(groups0)
    peak_pulse = np.maximum(pulse1, pulse2)
    print(f"     Peak pulse summary: min={np.min(peak_pulse)}, max={np.max(peak_pulse)}, mean={np.mean(peak_pulse):.2f}", flush=True)
    print("[12/14] Finish", flush=True)
    return peak_pulse


# In[ ]:


# ------------------------------------------------------------
# 模块 13: 计算评分并筛选 FP / Trash（支持 SNR 和 Score 两种模式）
# ------------------------------------------------------------
def compute_scores_and_filter(c_peak_mz_all, cell_occur, peak_pulse, aver_intensity,
                              blank_aver, blank_std, blank_rsd, rsd, n_cell_pulse,
                              filter_mode, filter_threshold):
    print(f"[13/14] MPCell scoring and cell-specific peak identification (cell pulse={n_cell_pulse})", flush=True)
    print(f"     Filter mode: {filter_mode}, filter threshold: {filter_threshold}", flush=True)
    cell_occur_ratio = cell_occur / n_cell_pulse
    base_mask = (cell_occur > 0.05 * n_cell_pulse)  # 细胞出现比例基本条件

    # 信噪比 sn
    sn = np.zeros_like(aver_intensity)
    mask_nonzero_std = (blank_std != 0)
    sn[mask_nonzero_std] = (aver_intensity[mask_nonzero_std] - blank_aver[mask_nonzero_std]) / blank_std[mask_nonzero_std]
    mask_std0_blank_non0 = (~mask_nonzero_std) & (blank_aver != 0)
    sn[mask_std0_blank_non0] = (aver_intensity[mask_std0_blank_non0] - blank_aver[mask_std0_blank_non0]) / blank_aver[mask_std0_blank_non0]
    mask_std0_blank0 = (~mask_nonzero_std) & (blank_aver == 0) & (aver_intensity != 0)
    sn[mask_std0_blank0] = 999

    pulse_ratio = np.full_like(peak_pulse, 999.0, dtype=float)
    mask_occur_nonzero = (cell_occur != 0)
    np.divide(peak_pulse, cell_occur, out=pulse_ratio, where=mask_occur_nonzero)

    rsd_subs = rsd - blank_rsd

    # 计算 total_score
    score_sn = np.zeros_like(sn)
    idx_high = np.where(sn >= 8)[0]
    score_sn[idx_high] = 10
    idx_mid = np.where((sn >= 5) & (sn < 8))[0]
    score_sn[idx_mid] = 10 * (sn[idx_mid] - 5) / 3

    score_rsd = np.zeros_like(rsd_subs)
    idx_high = np.where(rsd_subs >= 0.2)[0]
    score_rsd[idx_high] = 20
    idx_mid = np.where((rsd_subs >= 0.1) & (rsd_subs < 0.2))[0]
    score_rsd[idx_mid] = 10 * (rsd_subs[idx_mid] - 0.1) / 0.1

    score_pulse = np.zeros_like(pulse_ratio)
    idx_low = np.where(pulse_ratio <= 1.5)[0]
    score_pulse[idx_low] = 20
    idx_mid = np.where((pulse_ratio > 1.5) & (pulse_ratio <= 6))[0]
    score_pulse[idx_mid] = 10 * (6 - pulse_ratio[idx_mid]) / 4.5

    score_occur = np.zeros_like(cell_occur_ratio)
    idx_high = np.where(cell_occur_ratio >= 0.8)[0]
    score_occur[idx_high] = 20
    idx_mid = np.where((cell_occur_ratio >= 0.5) & (cell_occur_ratio < 0.8))[0]
    score_occur[idx_mid] = 10 * (cell_occur_ratio[idx_mid] - 0.5) / 0.3

    total_score = score_sn + 0.5 * (score_rsd + score_pulse + score_occur)
    #total_score=np.zeros_like(score_sn)
    #for idx,each in enumerate(total_score):
        #total_score[idx]=np.max([score_sn[idx],score_pulse[idx],score_occur[idx],score_rsd[idx]])

    # 根据模式筛选
    if filter_mode == 'snr':
        mask_fp = (sn >= filter_threshold) & base_mask
    else:  # 'score'
        mask_fp = (total_score >= filter_threshold) & base_mask

    fp_mz = c_peak_mz_all[mask_fp]
    trash_mz = c_peak_mz_all[~mask_fp]

    fp_props = {
        'mz': fp_mz,
        'pulse_ratio': pulse_ratio[mask_fp],
        'rsd': rsd_subs[mask_fp],
        'blank_rsd': blank_rsd[mask_fp],
        'cell_aver_intens': aver_intensity[mask_fp],
        'blank_aver_intens': blank_aver[mask_fp],
        'sn': sn[mask_fp],
        'cell_occur': cell_occur[mask_fp],
        'score': total_score[mask_fp]
    }
    trash_props = {
        'mz': trash_mz,
        'pulse_ratio': pulse_ratio[~mask_fp],
        'rsd': rsd_subs[~mask_fp],
        'blank_rsd': blank_rsd[~mask_fp],
        'cell_aver_intens': aver_intensity[~mask_fp],
        'blank_aver_intens': blank_aver[~mask_fp],
        'sn': sn[~mask_fp],
        'cell_occur': cell_occur[~mask_fp],
        'score': total_score[~mask_fp]
    }
    print(f"     Result: Positive = {len(fp_mz)} , Trash = {len(trash_mz)} ", flush=True)
    print("[13/14] Finish", flush=True)
    return fp_props, trash_props, sn, pulse_ratio, rsd_subs, total_score


# In[ ]:


# ------------------------------------------------------------
# 模块 14: 绘制 EIC 图（支持用户指定的时间范围）
# ------------------------------------------------------------
def plot_eics(path_p, RT_sub, eic_matrix_p, c_peak_mz_all, marker, ppm,
              fp_mz, fp_sn, trash_mz, trash_score, plot_option, rt_min_plot, rt_max_plot):
    print(f"[14/14] Draw EIC plot (plot option: {plot_option}, RT range: [{rt_min_plot}, {rt_max_plot}])", flush=True)
    if plot_option == 'none':
        print("     Skip plotting", flush=True)
        print("[14/14] Finish", flush=True)
        return

    # 找到 marker 的 EIC
    marker_idx = None
    for idx, mz_val in enumerate(c_peak_mz_all):
        if abs(mz_val - marker) < marker * ppm:
            marker_idx = idx
            break
    if marker_idx is None:
        print("     No marker EIC detected, skip", flush=True)
        print("[14/14] Finish", flush=True)
        return
    eic_marker = eic_matrix_p[marker_idx]

    # 绘制 FP
    if plot_option in ('fp', 'all'):
        path_fp = os.path.join(path_p, 'plot_fp')
        os.makedirs(path_fp, exist_ok=True)
        cnt = 0
        for idx, mz_val in enumerate(fp_mz):
            pos = np.where(np.abs(c_peak_mz_all - mz_val) < 1e-6)[0]
            if len(pos) == 0:
                continue
            peak_idx = pos[0]
            intensity = eic_matrix_p[peak_idx]
            plt.figure(figsize=(10,4))
            plt.subplot(2,1,1)
            plt.plot(RT_sub, eic_marker, label='marker')
            plt.xlim(rt_min_plot, rt_max_plot)
            plt.xlabel('RT/s')
            plt.ylabel('Int.')
            plt.legend()
            plt.subplot(2,1,2)
            plt.plot(RT_sub, intensity, label=f'm/z {mz_val:.4f}')
            plt.xlim(rt_min_plot, rt_max_plot)
            plt.xlabel('RT/s')
            plt.ylabel('Int.')
            plt.legend()
            plt.savefig(os.path.join(path_fp, f'snr{fp_sn[idx]:.2f} mz{mz_val:.4f}.png'))
            plt.close()
            cnt += 1
        print(f"     Drawn {cnt} EIC plots of positive peaks", flush=True)

    # 绘制 Trash
    if plot_option in ('trash', 'all'):
        path_trash = os.path.join(path_p, 'plot_trash')
        os.makedirs(path_trash, exist_ok=True)
        cnt = 0
        for idx, mz_val in enumerate(trash_mz):
            pos = np.where(np.abs(c_peak_mz_all - mz_val) < 1e-6)[0]
            if len(pos) == 0:
                continue
            peak_idx = pos[0]
            intensity = eic_matrix_p[peak_idx]
            plt.figure(figsize=(10,4))
            plt.subplot(2,1,1)
            plt.plot(RT_sub, eic_marker, label='marker')
            plt.xlim(rt_min_plot, rt_max_plot)
            plt.xlabel('RT/s')
            plt.ylabel('Int.')
            plt.legend()
            plt.subplot(2,1,2)
            plt.plot(RT_sub, intensity, label=f'm/z {mz_val:.4f}')
            plt.xlim(rt_min_plot, rt_max_plot)
            plt.xlabel('RT/s')
            plt.ylabel('Int.')
            plt.legend()
            plt.savefig(os.path.join(path_trash, f'score{trash_score[idx]:.2f} mz{mz_val:.4f}.png'))
            plt.close()
            cnt += 1
        print(f"     Drawn {cnt} EIC plots of trash", flush=True)

    print("[14/14] Finish", flush=True)


# In[ ]:


# ------------------------------------------------------------
# 主整合函数（需要接收绘图时间范围参数）
# ------------------------------------------------------------
def process_single_cell_data(params):
    print("\n========== Starting MPCell ==========", flush=True)
    # 打印参数摘要
    print("Parameter setting:", flush=True)
    print(f"  peak_intensity_thres = {params['peak_intensity_thres']}", flush=True)
    print(f"  ppm = {params['ppm']}", flush=True)
    print(f"  rt_range = [{params['rt_range_min']}, {params['rt_range_max']}]", flush=True)
    print(f"  cell_spec = {params['cell_spec']}, cell_spec_lim = {params['cell_spec_lim']}", flush=True)
    print(f"  marker = {params['marker']}, intensity_thres = {params['intensity_thres']}", flush=True)
    print(f"  bg_marker = {params['bg_marker']}", flush=True)
    print(f"  bg_intensity_max = {params['bg_intensity_max']}", flush=True)
    print(f"  input_file = {params['input_file']}", flush=True)
    print(f"  output_dir = {params['output_dir']}", flush=True)
    print(f"  plot_option = {params.get('plot_option', 'none')}", flush=True)
    print(f"  plot_rt_min = {params.get('plot_rt_min', 5)}", flush=True)
    print(f"  plot_rt_max = {params.get('plot_rt_max', 10)}", flush=True)
    print(f"  filter_mode = {params.get('filter_mode', 'score')}, filter_threshold = {params.get('filter_threshold', 10)}", flush=True)

    # 提取参数
    peak_intensity_thres = params['peak_intensity_thres']
    ppm = params['ppm']
    rt_range_min = params['rt_range_min']
    rt_range_max = params['rt_range_max']
    cell_spec = params['cell_spec']
    cell_spec_lim = params['cell_spec_lim']
    marker = params['marker']
    intensity_thres = params['intensity_thres']
    bg_marker = np.array(params['bg_marker'])
    bg_intensity_max = np.array(params['bg_intensity_max'])
    input_file = params['input_file']
    output_dir = params['output_dir']
    plot_option = params.get('plot_option', 'none')
    plot_rt_min = params.get('plot_rt_min', 5.0)
    plot_rt_max = params.get('plot_rt_max', 10.0)
    filter_mode = params.get('filter_mode', 'score')
    filter_threshold = params.get('filter_threshold', 10.0)

    os.makedirs(output_dir, exist_ok=True)

    # 1. 读取 mzML
    rt_raw, rt_to_ions = read_mzml_and_extract_ions(input_file, peak_intensity_thres, rt_range_max)
    np.save(os.path.join(output_dir, 'rt_to_ions.npy'), rt_to_ions)

    # 2. 筛选 RT 范围
    rt, index = filter_rt_range(rt_raw, rt_range_min, rt_range_max)
    np.save(os.path.join(output_dir, 'rt.npy'), rt)
    np.save(os.path.join(output_dir, 'index.npy'), index)

    # 3. 提取 marker 强度
    marker_intensity = extract_marker_intensity(rt, rt_to_ions, marker, ppm)
    df_marker = pd.DataFrame({'rt': rt, 'marker': marker_intensity})
    df_marker.to_excel(os.path.join(output_dir, 'marker_intensity.xlsx'), index=False)

    # 4. 检测细胞脉冲
    cell_pulse, cell_marker_vals, cell_rt, cell_index, cell_intensity = detect_cell_pulses(
        rt, index, marker_intensity, intensity_thres, cell_spec)
    if len(cell_pulse) == 0:
        print("No cell pulses detected, terminated", flush=True)
        return
    np.save(os.path.join(output_dir, 'cell_pulse.npy'), np.array(cell_pulse))
    with open(os.path.join(output_dir, 'cell_select_infor.txt'), 'w') as f:
        f.write(f'cell_num: {len(cell_pulse)}\n')
        f.write(f'cell_marker: {marker}\n')
        f.write(f'intensity_thres: {intensity_thres}\n')
    if len(cell_pulse) > 0:
        df_cell = pd.DataFrame({'cell-rt': [float(x[0]) for x in cell_pulse],
                                'marker': [np.max(x) for x in cell_marker_vals]})
        df_cell.to_excel(os.path.join(output_dir, 'cell_marker.xlsx'), index=False)

    # 5. 特征峰聚类
    cell_pulse_for_cluster = cell_pulse[:100] if len(cell_pulse) >= 100 else cell_pulse
    rt_list_for_cluster = sorted(set([rt0 for pulse in cell_pulse_for_cluster for rt0 in pulse]))
    mz_tar = cluster_features(rt_list_for_cluster, rt_to_ions, ppm, cell_pulse_for_cluster)
    np.save(os.path.join(output_dir, 'mz_tar.npy'), mz_tar)

    # 6. 构建 EIC 矩阵
    eic_matrix = extract_eic_matrix(rt, rt_to_ions, mz_tar, ppm)
    np.save(os.path.join(output_dir, 'c_peak_eic_intensity.npy'), eic_matrix)

    # 7. 提取细胞脉冲内的峰值信息
    c_peak = extract_cell_peak_max(rt, cell_pulse, eic_matrix, mz_tar)
    np.save(os.path.join(output_dir, 'c_peak.npy'), c_peak)

    # 8. 创建分组并提取子集
    path_p, cell_pulse_p, c_peak_p, eic_matrix_p, RT_sub = create_groups_and_subset(
        output_dir, rt, cell_pulse, c_peak, eic_matrix, bg_marker, bg_intensity_max, ppm, cell_spec_lim, marker)

    # 9. 平均强度和 RSD
    aver_intensity, rsd = compute_average_intensity_and_rsd(c_peak_p)
    np.save(os.path.join(path_p, 'c_peak_aver_intensity.npy'), aver_intensity)
    np.save(os.path.join(path_p, 'c_peak_rsd.npy'), rsd)

    # 10. 背景统计
    blank_aver, blank_std, blank_rsd, _ = compute_background_stats(
        path_p, bg_marker, bg_intensity_max, ppm, mz_tar, eic_matrix_p, cell_pulse_p, RT_sub)
    np.save(os.path.join(path_p, 'c_peak_blank_aver_intensity.npy'), blank_aver)
    np.save(os.path.join(path_p, 'c_peak_blank_std.npy'), blank_std)
    np.save(os.path.join(path_p, 'c_peak_blank_rsd.npy'), blank_rsd)

    # 11. 细胞出现次数
    cell_occur = compute_cell_occurrence(eic_matrix_p, cell_pulse_p, blank_aver, blank_std, RT_sub, cell_spec_lim)
    np.save(os.path.join(path_p, 'cell_occur.npy'), cell_occur)

    # 12. 峰脉冲数
    peak_pulse = compute_peak_pulse_count(eic_matrix_p, blank_aver, blank_std)
    np.save(os.path.join(path_p, 'c_peak_pulse.npy'), peak_pulse)

    # 13. 评分及筛选
    fp_props, trash_props, sn, pulse_ratio, rsd_subs, total_score = compute_scores_and_filter(
        mz_tar, cell_occur, peak_pulse, aver_intensity, blank_aver, blank_std, blank_rsd, rsd, len(cell_pulse_p),
        filter_mode, filter_threshold)
    np.save(os.path.join(path_p, 'c_peak_sn.npy'), sn)
    np.save(os.path.join(path_p, 'c_peak_pulse_ratio.npy'), pulse_ratio)
    np.save(os.path.join(path_p, 'c_peak_rsd_subs.npy'), rsd_subs)

    # 保存 Excel 等（与之前相同）
    df_fp = pd.DataFrame(fp_props)
    df_fp.to_excel(os.path.join(path_p, 'c_peak_fp.xlsx'), index=False)
    df_trash = pd.DataFrame(trash_props)
    df_trash.to_excel(os.path.join(path_p, 'c_peak_trash.xlsx'), index=False)
    df_all = pd.DataFrame({
        'm/z': mz_tar,
        'pulse_ratio': pulse_ratio,
        'rsd': rsd_subs,
        'blank_rsd': blank_rsd,
        'cell_aver_intens': aver_intensity,
        'blank_aver_intens': blank_aver,
        'sn': sn,
        'cell_occur': cell_occur
    })
    df_all.to_excel(os.path.join(path_p, 'c_peak_infor.xlsx'), index=False)

    np.save(os.path.join(path_p, 'c_peak_fp.npy'), fp_props['mz'])
    np.save(os.path.join(path_p, 'c_peak_trash.npy'), trash_props['mz'])

    # 14. 绘图（传入绘图时间范围）
    plot_eics(path_p, RT_sub, eic_matrix_p, mz_tar, marker, ppm,
              fp_props['mz'], fp_props['sn'], trash_props['mz'], trash_props['score'],
              plot_option, plot_rt_min, plot_rt_max)

    # 最终汇总
    print("\n========== END ==========", flush=True)
    print(f"cell pulse: {len(cell_pulse)}", flush=True)
    print(f"mz_list length: {len(mz_tar)}", flush=True)
    print(f"cell-specific peak: {len(fp_props['mz'])}", flush=True)
    print(f"trash peak: {len(trash_props['mz'])}", flush=True)
    print(f"All results saved in: {output_dir}", flush=True)
    print(f"All subgroup results saved in: {path_p}", flush=True)

# ------------------------------------------------------------
# 主入口
# ------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Error：file path not specified", file=sys.stderr, flush=True)
        sys.exit(1)
    param_file = sys.argv[1]
    try:
        with open(param_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Failed to read file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)
    try:
        process_single_cell_data(data)
    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr, flush=True)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()


# In[ ]:




