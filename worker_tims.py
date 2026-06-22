#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
    # 读取第一个谱图获取 CCS 参考值
    reader = mzml.read(input_file)
    ccs_spec = []
    for spectrum in reader:
        rt_ini = float(spectrum['scanList']['scan'][0]['scan start time'])
        break
    reader = mzml.read(input_file)
    for spectrum in reader:
        rt0 = float(spectrum['scanList']['scan'][0]['scan start time'])
        ccs0 = float(spectrum['scanList']['scan'][0]['inverse reduced ion mobility'])
        if rt0 == rt_ini:
            ccs_spec.append(ccs0)
        else:
            break

    rt_to_ions = {}
    reader = mzml.read(input_file)
    rt = []
    for spectrum in reader:
        rt0 = float(spectrum['scanList']['scan'][0]['scan start time'])
        ccs0 = float(spectrum['scanList']['scan'][0]['inverse reduced ion mobility'])
        if ccs0 == ccs_spec[0]:
            rt.append(rt0)
        mz_array = np.array(spectrum['m/z array'])
        intensity_array = np.array(spectrum['intensity array'])
        mask = intensity_array > peak_intensity_thres
        if not np.any(mask):
            continue
        mz_array = mz_array[mask]
        intensity_array = intensity_array[mask]
        ccs_array = np.full_like(mz_array, ccs0)
        if rt0 in rt_to_ions:
            prev_mz, prev_ccs, prev_int = rt_to_ions[rt0]
            new_mz = np.concatenate((prev_mz, mz_array))
            new_ccs = np.concatenate((prev_ccs, ccs_array))
            new_int = np.concatenate((prev_int, intensity_array))
            order = np.argsort(new_mz)
            rt_to_ions[rt0] = (new_mz[order], new_ccs[order], new_int[order])
        else:
            order = np.argsort(mz_array)
            rt_to_ions[rt0] = (mz_array[order], ccs_array[order], intensity_array[order])
        if rt0 > rt_range_max:
            break
    print(f"     Read {len(rt)} scans (MS1 spec)，effective spec: {len(rt_to_ions)}", flush=True)
    print("[1/14] Finish", flush=True)
    return rt, rt_to_ions, ccs_spec

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
def extract_marker_intensity(rt, rt_to_ions, marker, mzppm, ccsppm):
    print(f"[3/14] Extract marker m/z = {marker[1]} (mzppm={mzppm}), CCS={marker[0]} (ccsppm={ccsppm}) EIC", flush=True)
    marker_mz = marker[1]
    marker_ccs = marker[0]
    marker_intensity = np.zeros(len(rt))
    detected_count = 0
    for idx, rt0 in enumerate(rt):
        mz_arr, ccs_arr, intensity_arr = rt_to_ions[rt0]
        mz_tol = marker_mz * mzppm
        ccs_tol = marker_ccs * ccsppm
        mask = (abs(mz_arr - marker_mz) < mz_tol) & (abs(ccs_arr - marker_ccs) < ccs_tol)
        if np.any(mask):
            marker_intensity[idx] = np.sum(intensity_arr[mask])
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


# In[ ]:


# ------------------------------------------------------------
# 模块 5: TIMS 数据特征峰聚类（基于 CCS 维度，区域增长 + 堆排序）
# ------------------------------------------------------------
def cluster_features_tims(ccs_spec, rt_to_ions, cell_pulse, mzppm, eim_span=27):
    """
    参数:
        path: 输出目录（保存中间结果）
        ccs_spec: 一维数组，所有可能的 CCS 值
        rt_to_ions: 字典，键为 rt0，值为 (mz_arr, ccs_arr, intensity_arr)
        cell_pulse: 细胞脉冲列表，每个脉冲包含一系列 RT 点
        mzppm: m/z 容差（ppm）
        eim_span: 精筛时的窗口宽度（默认 27）
    返回:
        peak_tar_filter_arr: 筛选后的特征峰数组，每个元素为 (CCS, mz)
    """
    import bisect
    import heapq
    import numpy as np
    from statsmodels.nonparametric.smoothers_lowess import lowess
    from scipy.signal import find_peaks

    print("[5/14] TIMS peak grouping", flush=True)

    # ---------- 步骤1: 构建每个 CCS 值对应的离子列表 ----------
    peak_arr = np.zeros(len(ccs_spec), dtype=object)
    for idx, each_ccs in enumerate(ccs_spec):
        peak_arr_each = []
        for pulse in cell_pulse:
            for rt0 in pulse:
                mz_arr, ccs_arr, intensity_arr = rt_to_ions[rt0]
                mask = np.where(ccs_arr == each_ccs)[0]
                mz_arr_each = mz_arr[mask]
                intensity_arr_each = intensity_arr[mask]
                for x in range(len(mz_arr_each)):
                    peak_arr_each.append((mz_arr_each[x], intensity_arr_each[x]))
        peak_arr[idx] = np.array(peak_arr_each)
    print("     Peak list constructed", flush=True)

    # ---------- 步骤2: 对每个 CCS 值进行 m/z 聚类（迭代聚类）----------
    cell_aggre = np.zeros(len(peak_arr), dtype=object)
    for idx0, mz_list in enumerate(peak_arr):
        if idx0 % 100 == 0:
            print(f"     Schedule: {idx0}/{len(peak_arr)}", flush=True)
        if not np.any(mz_list):
            continue
        mz_list = mz_list.tolist()
        mz_list.sort(key=lambda x: x[1], reverse=True)
        mz_list_arr = np.array(mz_list)

        cycle = 0
        while True:
            clusters = []
            ref_mzs = []
            for (mz, intensity) in mz_list_arr:
                pos = bisect.bisect_left(ref_mzs, mz)
                found = False
                for candidate_idx in (pos, pos - 1):
                    if 0 <= candidate_idx < len(clusters):
                        ref = clusters[candidate_idx]['ref_mz']
                        if abs(mz - ref) <= ref * mzppm:
                            clusters[candidate_idx]['points'].append((mz, intensity))
                            found = True
                            break
                if not found:
                    new_cluster = {'ref_mz': mz, 'points': [(mz, intensity)]}
                    clusters.insert(pos, new_cluster)
                    ref_mzs.insert(pos, mz)

            percent = 0.01 if cycle == 0 else 0
            mz_list_final = []
            for cluster in clusters:
                points = np.array(cluster['points'])
                if len(points) > percent * len(cell_pulse):
                    plus_sum = np.sum(points[:, 0] * points[:, 1])
                    intensity_sum = np.sum(points[:, 1])
                    mz_tar = plus_sum / intensity_sum
                    mz_list_final.append((mz_tar, intensity_sum))
            mz_list_final.sort(key=lambda x: x[1], reverse=True)
            mz_list_arr = np.array(mz_list_final, dtype=object)
            cycle += 1
            if len(mz_list) == len(mz_list_arr):
                break
            mz_list = mz_list_arr.tolist()
        cell_aggre[idx0] = mz_list_arr
    print("     CCS-based grouping finished", flush=True)

    # ---------- 步骤3: 构建 cell_super_aggre（添加 CCS 信息）----------
    cell_super_aggre = []
    for idx0, each in enumerate(cell_aggre):
        temp = []
        ccs0 = ccs_spec[idx0]
        if np.any(each):
            for peak in each:
                temp.append([ccs0] + peak.tolist())
            cell_super_aggre.append(temp)
        else:
            cell_super_aggre.append([])
    cell_super_aggre = np.array(cell_super_aggre, dtype=object)
    print("     cell_super_aggre constructed", flush=True)

    # ---------- 步骤4: 堆排序区域增长（提取 EIM）----------
    heap = []
    total_ions = 0
    for scan_idx, scan_ions in enumerate(cell_super_aggre):
        for ion_idx, (ccs, mz, ints) in enumerate(scan_ions):
            heapq.heappush(heap, (-ints, scan_idx, ion_idx))
            total_ions += 1

    peak_tar = []
    processed_regions = 0
    while total_ions > 0:
        if not heap:
            break

        if processed_regions > 0 and processed_regions % 50000 == 0:
            print('     Rebuild heap...', end='', flush=True)
            cell_super_aggre_no0 = []
            for scan_idx, scan_ions in enumerate(cell_super_aggre):
                temp = []
                for ion_idx, (ccs, mz, ints) in enumerate(scan_ions):
                    if ints != 0:
                        temp.append(cell_super_aggre[scan_idx][ion_idx])
                cell_super_aggre_no0.append(temp)
            cell_super_aggre = cell_super_aggre_no0
            new_heap = []
            for scan_idx, scan_ions in enumerate(cell_super_aggre):
                for ion_idx, (ccs, mz, ints) in enumerate(scan_ions):
                    heapq.heappush(new_heap, (-ints, scan_idx, ion_idx))
            heap = new_heap
            print('ok', flush=True)
            if not heap:
                break

        processed_regions += 1
        neg_ints, scan_idx, ion_idx = heapq.heappop(heap)
        if scan_idx >= len(cell_super_aggre) or ion_idx >= len(cell_super_aggre[scan_idx]):
            continue
        if cell_super_aggre[scan_idx][ion_idx][-1] == 0 or cell_super_aggre[scan_idx][ion_idx][-1] != -neg_ints:
            continue

        ccs0, mz0, int0 = cell_super_aggre[scan_idx][ion_idx]
        mz_tol_abs = mz0 * mzppm
        eim = []

        # 向左扩展
        for i in range(scan_idx, -1, -1):
            scan_ions = cell_super_aggre[i]
            if not scan_ions:
                break
            found = False
            for j, (ccs, mz, ints) in enumerate(scan_ions):
                if ints == 0:
                    continue
                if abs(mz - mz0) <= mz_tol_abs:
                    eim.append((ccs, mz, ints))
                    scan_ions[j] = (ccs, mz, 0)
                    total_ions -= 1
                    found = True
                    break
            if not found:
                break

        # 向右扩展
        for i in range(scan_idx + 1, len(cell_super_aggre)):
            scan_ions = cell_super_aggre[i]
            if not scan_ions:
                break
            found = False
            for j, (ccs, mz, ints) in enumerate(scan_ions):
                if ints == 0:
                    continue
                if abs(mz - mz0) <= mz_tol_abs:
                    eim.append((ccs, mz, ints))
                    scan_ions[j] = (ccs, mz, 0)
                    total_ions -= 1
                    found = True
                    break
            if not found:
                break

        if len(eim) >= 15:
            peak_tar.append(eim)
            print(f'     Peak targets: {len(peak_tar)}, total_ions: {total_ions}', flush=True)

    peak_tar = np.array(peak_tar, dtype=object)
    print("     {} candidate EIM extracted".format(len(peak_tar)), flush=True)

    # ---------- 步骤5: 精筛（基于平滑峰形和 SD）使用 eim_span 作为窗口宽度----------
    window_width = eim_span   # 使用传入的 eim_span，默认 27
    window = 10
    peak_tar_filter = []

    for idx0, eim_part in enumerate(peak_tar):
        if idx0 % 100 == 0:
            print(f"     Schedule for precise filtration: {idx0}/{len(peak_tar)}", flush=True)
        eim_part.sort(key=lambda x: x[0], reverse=True)
        eim_part_arr = np.array(eim_part)

        eim = np.zeros(len(ccs_spec))
        for (ccs, mz, intensity) in eim_part_arr:
            for i in range(len(eim)):
                if ccs == ccs_spec[i]:
                    eim[i] = intensity
                    break
        raw_intensity = eim

        # 平滑
        frac = window / len(raw_intensity)
        smoothed = lowess(raw_intensity, np.arange(len(raw_intensity)), frac=frac, it=0, delta=0.0, return_sorted=False)

        # 找顶点
        mask = np.where(eim != 0)[0]
        if len(mask) == 0:
            continue
        start = int(mask[0])
        peaks, _ = find_peaks(smoothed[mask])
        if len(peaks) == 0:
            continue
        apex = start + peaks[np.argmax(smoothed[start + peaks])]

        # 局部窗口
        half = window_width // 2
        left = max(0, apex - half)
        right = min(len(raw_intensity) - 1, apex + half + 1)
        raw_win = raw_intensity[left:right]
        smooth_win = smoothed[left:right]

        # 计算 peak_sd
        diff = raw_win - smooth_win
        sd = np.std(diff, ddof=1)
        smooth_range = np.max(smooth_win) - np.min(smooth_win)
        if smooth_range == 0:
            continue
        peak_sd = sd / smooth_range

        if peak_sd < 0.35:
            plus_sum = 0
            intensity_sum = 0
            for (ccs, mz, intensity) in eim_part_arr:
                plus_sum += mz * intensity
                intensity_sum += intensity
            if intensity_sum != 0:
                mz_tar = plus_sum / intensity_sum
                peak_tar_filter.append((ccs_spec[apex], mz_tar))
                print(f"     Filtered peak: {len(peak_tar_filter)}", flush=True)

    peak_tar_filter_arr = np.array(peak_tar_filter, dtype=object)
    print(f"     Precise filtration finished, {len(peak_tar_filter)} peak obtained (CCS, mz)", flush=True)
    return peak_tar_filter_arr


# In[ ]:


# ------------------------------------------------------------
# 模块 6: 构建 c_peak EIC 强度矩阵（支持 CCS 容差）
# ------------------------------------------------------------
def extract_eic_matrix_tims(rt, rt_to_ions, peak_tar_filter, mzppm, ccsppm):
    """
    参数:
        rt: 一维保留时间数组
        rt_to_ions: 字典 {rt0: (mz_arr, ccs_arr, intensity_arr)}
        peak_tar_filter: numpy 数组，形状 (n_peaks, 2)，每行为 [ccs, mz]
        mzppm: m/z 容差 (ppm)
        ccsppm: CCS 容差 (ppm)
    返回:
        time_intensity_orig: 形状 (n_peaks, len(rt)) 的强度矩阵
    """
    print(f"[6/14] Constrction of EIC intensity matrix (peak length={len(peak_tar_filter)}, time point={len(rt)})", flush=True)
    
    # 按 m/z 排序，提高二分查找效率
    sorted_idx = np.argsort(peak_tar_filter[:, 1])
    peak_tar_filter_sorted = peak_tar_filter[sorted_idx]
    orig_order = np.argsort(sorted_idx)
    
    n_times = len(rt)
    n_peaks = len(peak_tar_filter)
    time_intensity = np.zeros((n_peaks, n_times), dtype=np.float32)
    
    for t_idx, rt0 in enumerate(rt):
        if t_idx % max(1, n_times // 10) == 0:
            print(f"     Processed time points {t_idx}/{n_times}", flush=True)
        if rt0 not in rt_to_ions:
            continue
        mz_arr, ccs_arr, int_arr = rt_to_ions[rt0]
        # 对每个峰（按 m/z 排序）查找匹配离子
        for peak_idx in range(n_peaks):
            mz_target = peak_tar_filter_sorted[peak_idx, 1]
            ccs_target = peak_tar_filter_sorted[peak_idx, 0]
            tol = mz_target * mzppm
            left = bisect.bisect_left(mz_arr, mz_target - tol)
            right = bisect.bisect_right(mz_arr, mz_target + tol)
            if left >= right:
                continue
            # 在候选区间内筛选 CCS
            ccs_candidate = ccs_arr[left:right]
            int_candidate = int_arr[left:right]
            mask_ccs = np.abs(ccs_candidate - ccs_target) < ccsppm
            if np.any(mask_ccs):
                time_intensity[peak_idx, t_idx] = np.sum(int_candidate[mask_ccs])
    
    # 恢复原始峰顺序
    time_intensity_orig = time_intensity[orig_order, :]
    print(f"     EIC shape: {time_intensity_orig.shape}", flush=True)
    print("[6/14] Finish", flush=True)
    return time_intensity_orig


# In[ ]:


# ------------------------------------------------------------
# 模块 7: 提取每个细胞脉冲内峰的最大强度和对应 RT
# ------------------------------------------------------------
def extract_cell_peak_max_tims(rt, cell_pulse, eic_matrix, peak_tar_filter):
    """
    参数:
        rt: 一维保留时间数组
        cell_pulse: 细胞脉冲列表，每个元素为该脉冲包含的 RT 列表
        eic_matrix: 形状 (n_peaks, len(rt)) 的强度矩阵
        peak_tar_filter: numpy 数组，形状 (n_peaks, 2)，每行为 [ccs, mz]
    返回:
        c_peak_max_arr: object 数组，每个元素为 [ [ccs, mz], rt_max_list, intensity_max_list ]
    """
    print(f"[7/14] Max intensity extraction in cell pulses ({len(cell_pulse)} cell pulses, {len(peak_tar_filter)} peaks)", flush=True)
    
    rt_to_index = {rt0: idx for idx, rt0 in enumerate(rt)}
    n_peaks = len(peak_tar_filter)
    c_peak_max = []
    
    for each in range(n_peaks):
        if each % max(1, n_peaks // 10) == 0:
            print(f"     Schedule: {each}/{n_peaks}", flush=True)
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
        c_peak_max.append([peak_tar_filter[each], rt_max_list, intensity_max_list])
    
    c_peak_max_arr = np.array(c_peak_max, dtype=object)
    print("[7/14] Finish", flush=True)
    return c_peak_max_arr


# In[ ]:


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
    print("[8/14] 完成", flush=True)
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
    print(f"     Finish, averaged intensity rang: [{np.min(aver_intensity):.2e}, {np.max(aver_intensity):.2e}], averaged RSD: {np.mean(rsd):.4f}", flush=True)
    print("[9/14] Finish", flush=True)
    return aver_intensity, rsd

# ------------------------------------------------------------
# 模块 10: 计算背景（非细胞区域）的平均强度、标准差、RSD（TIMS 版本）
# ------------------------------------------------------------
def compute_background_stats_tims(path_p, bg_marker, bg_intensity_max, mzppm, ccsppm,
                                  peak_tar_filter, eic_matrix_p, cell_pulse_p, RT_sub):
    """
    参数:
        path_p: 分组子目录路径
        bg_marker: 背景标记列表，每个元素为 (ccs, mz)
        bg_intensity_max: 背景最大强度列表
        mzppm: m/z 容差（ppm）
        ccsppm: CCS 容差（ppm）
        peak_tar_filter: 特征峰数组，形状 (n_peaks, 2)，每行 [ccs, mz]
        eic_matrix_p: 当前分组的 EIC 矩阵 (n_peaks, n_timepoints)
        cell_pulse_p: 当前分组的细胞脉冲列表
        RT_sub: 当前分组的时间点数组
    返回:
        blank_aver, blank_std, blank_rsd, c_peak_blank_index
    """
    print(f"[10/14] Peak information in background (bg_marker={bg_marker}, bg_intensity_max={bg_intensity_max})", flush=True)
    
    # 找出所有背景时间点（不在任何细胞脉冲内）
    cell_pulse_rt = set()
    for pulse in cell_pulse_p:
        cell_pulse_rt.update(pulse)
    cell_pulse_rt = sorted(cell_pulse_rt)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(RT_sub)}
    cell_pulse_idx = [rt_to_index[rt0] for rt0 in cell_pulse_rt if rt0 in rt_to_index]
    bg_nocell_idx = [idx for idx in range(len(RT_sub)) if idx not in cell_pulse_idx]
    print(f"     no-cell BG time point: {len(bg_nocell_idx)} / {len(RT_sub)}", flush=True)

    # 根据每个背景标记筛选符合强度条件的背景时间点
    c_peak_blank_index = None
    for count, (bg_mz_ccs, bg_max) in enumerate(zip(bg_marker, bg_intensity_max)):
        # bg_marker 每个元素是 (ccs, mz)
        bg_ccs, bg_mz = bg_mz_ccs
        # 寻找与背景标记最接近的峰索引（同时考虑 m/z 和 CCS）
        mask = (np.abs(peak_tar_filter[:, 1] - bg_mz) < bg_mz * mzppm) & (np.abs(peak_tar_filter[:, 0] - bg_ccs) < ccsppm)
        peak_indices = np.where(mask)[0]
        if len(peak_indices) == 0:
            print(f"     No peak matches background marker {bg_mz_ccs}, skip", flush=True)
            continue
        # 取第一个匹配的峰
        eic_intensity_marker = eic_matrix_p[peak_indices[0]]
        bg_intensity_idx = np.where(eic_intensity_marker <= bg_max)[0]
        if count == 0:
            c_peak_blank_index = set(bg_intensity_idx)
        else:
            c_peak_blank_index.intersection_update(bg_intensity_idx)
        print(f"     BG marker {bg_mz_ccs} filtered BG time point: {len(c_peak_blank_index)}", flush=True)

    if c_peak_blank_index is None:
        raise ValueError("No effective blank index，please check parameters of bg_marker")
    
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
        # 动态确定最小样本数，避免过小样本导致的不稳定
        min_samples = min(max(10, int(0.01 * len(c_peak_blank_index))), 50)
        if len(non_zero) > min_samples:
            sorted_vals = np.sort(non_zero)
            trim_len = int(0.1 * len(sorted_vals))
            trimmed = sorted_vals[trim_len:-trim_len] if len(sorted_vals) > 2 * trim_len else sorted_vals
        else:
            sorted_vals = np.sort(blank_intensity)
            trim_len = int(0.1 * len(sorted_vals))
            trimmed = sorted_vals[trim_len:-trim_len] if len(sorted_vals) > 2 * trim_len else sorted_vals
        mean_val = np.mean(trimmed)
        blank_aver[i] = mean_val
        blank_std[i] = np.std(trimmed, ddof=1)
        if mean_val != 0:
            blank_rsd[i] = blank_std[i] / mean_val
    
    print(f"     avergaed intensity range: [{np.min(blank_aver):.2e}, {np.max(blank_aver):.2e}]", flush=True)
    print("[10/14] 完成", flush=True)
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
    print("[11/14] 完成", flush=True)
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
# 模块 13: 计算评分并筛选 FP / Trash（TIMS 版本，支持 CCS 输出）
# ------------------------------------------------------------
def compute_scores_and_filter_tims(peak_tar_filter, cell_occur, peak_pulse, aver_intensity,
                                   blank_aver, blank_std, blank_rsd, rsd, n_cell_pulse,
                                   filter_mode, filter_threshold):
    """
    参数:
        peak_tar_filter: numpy 数组，形状 (n_peaks, 2)，每行 [ccs, mz]
        cell_occur: 细胞出现次数数组
        peak_pulse: 脉冲数数组
        aver_intensity: 细胞脉冲平均强度数组
        blank_aver: 背景平均强度数组
        blank_std: 背景标准差数组
        blank_rsd: 背景 RSD 数组
        rsd: 细胞脉冲 RSD 数组
        n_cell_pulse: 总细胞脉冲数
        filter_mode: 'snr' 或 'score'
        filter_threshold: 阈值
    返回:
        fp_props: dict, 包含 fp 的所有属性 (mz, ccs, pulse_ratio, rsd, blank_rsd,
                    cell_aver_intens, blank_aver_intens, sn, cell_occur, score)
        trash_props: 同理
        sn, pulse_ratio, rsd_subs, total_score
    """
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

    # 评分计算 (与用户提供的公式完全一致)
    score_sn = np.zeros_like(sn)
    idx_high = np.where(sn >= 8)[0]
    score_sn[idx_high] = 10
    idx_mid = np.where((sn >= 5) & (sn < 8))[0]
    score_sn[idx_mid] = 10 * (sn[idx_mid] - 5) / 3

    score_rsd = np.zeros_like(rsd_subs)
    idx_high = np.where(rsd_subs >= 0.2)[0]
    score_rsd[idx_high] = 10
    idx_mid = np.where((rsd_subs >= 0.1) & (rsd_subs < 0.2))[0]
    score_rsd[idx_mid] = 10 * (rsd_subs[idx_mid] - 0.1) / 0.1

    score_pulse = np.zeros_like(pulse_ratio)
    idx_low = np.where(pulse_ratio <= 1.5)[0]
    score_pulse[idx_low] = 10
    idx_mid = np.where((pulse_ratio > 1.5) & (pulse_ratio <= 6))[0]
    score_pulse[idx_mid] = 10 * (6 - pulse_ratio[idx_mid]) / 4.5

    score_occur = np.zeros_like(cell_occur_ratio)
    idx_high = np.where(cell_occur_ratio >= 0.8)[0]
    score_occur[idx_high] = 10
    idx_mid = np.where((cell_occur_ratio >= 0.5) & (cell_occur_ratio < 0.8))[0]
    score_occur[idx_mid] = 10 * (cell_occur_ratio[idx_mid] - 0.5) / 0.3

    total_score = score_sn + 0.5 * (score_rsd + score_pulse + score_occur)

    # 根据模式筛选
    if filter_mode == 'snr':
        mask_fp = (sn >= filter_threshold) & base_mask
    else:
        mask_fp = (total_score >= filter_threshold) & base_mask

    # 提取 FP 和 Trash 的 CCS + m/z
    fp_mz_ccs = peak_tar_filter[mask_fp]
    trash_mz_ccs = peak_tar_filter[~mask_fp]

    # 构建返回的字典
    fp_props = {
        'mz': fp_mz_ccs[:, 1],
        'ccs': fp_mz_ccs[:, 0],
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
        'mz': trash_mz_ccs[:, 1],
        'ccs': trash_mz_ccs[:, 0],
        'pulse_ratio': pulse_ratio[~mask_fp],
        'rsd': rsd_subs[~mask_fp],
        'blank_rsd': blank_rsd[~mask_fp],
        'cell_aver_intens': aver_intensity[~mask_fp],
        'blank_aver_intens': blank_aver[~mask_fp],
        'sn': sn[~mask_fp],
        'cell_occur': cell_occur[~mask_fp],
        'score': total_score[~mask_fp]
    }

    print(f"     Result: Positive = {len(fp_mz_ccs)}, Trash = {len(trash_mz_ccs)}", flush=True)
    print("[13/14] Finish", flush=True)
    return fp_props, trash_props, sn, pulse_ratio, rsd_subs, total_score


# In[ ]:


# ------------------------------------------------------------
# 模块 14: 绘制 EIC 图（TIMS 版本，显示 CCS 和 m/z）
# ------------------------------------------------------------
def plot_eics_tims(path_p, RT_sub, eic_matrix_p, peak_tar_filter, marker,
                   mzppm, ccsppm, fp_props, trash_props, plot_option,
                   rt_min_plot, rt_max_plot):
    """
    参数:
        path_p: 分组子目录路径
        RT_sub: 当前分组的时间点数组
        eic_matrix_p: EIC 矩阵 (n_peaks, n_timepoints)
        peak_tar_filter: 所有特征峰的 (CCS, mz) 数组，形状 (n_peaks, 2)
        marker: 用户设置的 marker (ccs, mz) 元组
        mzppm: m/z 容差 (ppm)
        ccsppm: CCS 容差 (ppm)
        fp_props: 包含 'ccs', 'mz', 'sn' 等字段的字典
        trash_props: 类似，包含 'ccs', 'mz', 'score'
        plot_option: 'none', 'fp', 'trash', 'all'
        rt_min_plot: 绘图 RT 下限
        rt_max_plot: 绘图 RT 上限
    """
    print(f"[14/14] 绘制 EIC 图 (选项: {plot_option}, RT 范围: [{rt_min_plot}, {rt_max_plot}])", flush=True)
    if plot_option == 'none':
        print("     Skip plotting", flush=True)
        print("[14/14] Finish", flush=True)
        return

    # 1. 找到 marker 对应的 EIC 强度序列
    marker_idx = None
    marker_ccs, marker_mz = marker
    for idx, (ccs, mz) in enumerate(peak_tar_filter):
        if (abs(mz - marker_mz) < marker_mz * mzppm) and (abs(ccs - marker_ccs) < ccsppm):
            marker_idx = idx
            break
    if marker_idx is None:
        print("     No marker EIC detected, skip", flush=True)
        print("[14/14] Finish", flush=True)
        return
    eic_marker = eic_matrix_p[marker_idx]

    # 2. 辅助函数：绘制单个 EIC 图
    def plot_single_eic(eic_intensity, title, filename, sn_or_score=None):
        plt.figure(figsize=(10, 4))
        # 子图1: marker
        plt.subplot(2, 1, 1)
        plt.plot(RT_sub, eic_marker, label='marker', color='blue')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Intensity')
        plt.legend()
        # 子图2: 当前峰
        plt.subplot(2, 1, 2)
        plt.plot(RT_sub, eic_intensity, label=title, color='red')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Intensity')
        plt.legend()
        # 保存图形
        plt.savefig(os.path.join(path_p, filename))
        plt.close()

    # 3. 绘制 FP
    if plot_option in ('fp', 'all'):
        path_fp = os.path.join(path_p, 'plot_fp')
        os.makedirs(path_fp, exist_ok=True)
        cnt = 0
        for i in range(len(fp_props['mz'])):
            ccs = fp_props['ccs'][i]
            mz = fp_props['mz'][i]
            sn = fp_props['sn'][i]
            # 找到该峰在 EIC 矩阵中的索引
            pos = np.where((np.abs(peak_tar_filter[:, 0] - ccs) < 1e-6) &
                           (np.abs(peak_tar_filter[:, 1] - mz) < 1e-6))[0]
            if len(pos) == 0:
                continue
            peak_idx = pos[0]
            intensity = eic_matrix_p[peak_idx]
            title = f'CCS={ccs:.4f}  m/z={mz:.4f}'
            filename = f'snr_{sn:.2f}_ccs{ccs:.4f}_mz{mz:.4f}.png'
            plot_single_eic(intensity, title, os.path.join('plot_fp', filename), sn)
            cnt += 1
        print(f"     Drawn {cnt} EIC plots of positive peaks", flush=True)

    # 4. 绘制 Trash
    if plot_option in ('trash', 'all'):
        path_trash = os.path.join(path_p, 'plot_trash')
        os.makedirs(path_trash, exist_ok=True)
        cnt = 0
        for i in range(len(trash_props['mz'])):
            ccs = trash_props['ccs'][i]
            mz = trash_props['mz'][i]
            score = trash_props['score'][i]
            pos = np.where((np.abs(peak_tar_filter[:, 0] - ccs) < 1e-6) &
                           (np.abs(peak_tar_filter[:, 1] - mz) < 1e-6))[0]
            if len(pos) == 0:
                continue
            peak_idx = pos[0]
            intensity = eic_matrix_p[peak_idx]
            title = f'CCS={ccs:.4f}  m/z={mz:.4f}'
            filename = f'score_{score:.2f}_ccs{ccs:.4f}_mz{mz:.4f}.png'
            plot_single_eic(intensity, title, os.path.join('plot_trash', filename), score)
            cnt += 1
        print(f"     Drawn {cnt} EIC plots of trash", flush=True)

    print("[14/14] Finish", flush=True)


# In[ ]:


# ------------------------------------------------------------
# 主整合函数（需要接收绘图时间范围参数）
# ------------------------------------------------------------
def process_single_cell_data(params):
    print("\n========== Starting MPCell ==========", flush=True)
   
    # 提取参数（tims 模式专用）
    peak_intensity_thres = params['peak_intensity_thres']
    mzppm = params.get('mzppm', 20e-6)      # 从 GUI 传入，默认 20e-6
    ccsppm = params.get('ccsppm', 0.05)     # 从 GUI 传入，默认 0.05
    rt_range_min = params['rt_range_min']
    rt_range_max = params['rt_range_max']
    cell_spec = params['cell_spec']
    cell_spec_lim = params['cell_spec_lim']
    marker = params['marker']                # 应为 (ccs, mz) 元组
    intensity_thres = params['intensity_thres']
    bg_marker = params['bg_marker']          # 列表，每个元素 (ccs, mz)
    bg_intensity_max = params['bg_intensity_max']
    input_file = params['input_file']
    output_dir = params['output_dir']
    plot_option = params.get('plot_option', 'none')
    plot_rt_min = params.get('plot_rt_min', 5.0)
    plot_rt_max = params.get('plot_rt_max', 10.0)
    filter_mode = params.get('filter_mode', 'score')
    filter_threshold = params.get('filter_threshold', 10.0)
    eim_span = params.get('eim_span', 27)
    
    # 打印参数摘要
    print("Parameter setting:", flush=True)
    print(f"  peak_intensity_thres = {peak_intensity_thres}", flush=True)
    print(f"  mzppm = {mzppm}", flush=True)
    print(f"  ccsppm = {ccsppm}", flush=True)
    print(f"  rt_range = [{rt_range_min}, {rt_range_max}]", flush=True)
    print(f"  cell_spec = {cell_spec}, cell_spec_lim = {cell_spec_lim}", flush=True)
    print(f"  marker (ccs, mz) = {marker}, intensity_thres = {intensity_thres}", flush=True)
    print(f"  bg_marker = {bg_marker}", flush=True)
    print(f"  bg_intensity_max = {bg_intensity_max}", flush=True)
    print(f"  input_file = {input_file}", flush=True)
    print(f"  output_dir = {output_dir}", flush=True)
    print(f"  plot_option = {plot_option}", flush=True)
    print(f"  plot_rt_min = {plot_rt_min}, plot_rt_max = {plot_rt_max}", flush=True)
    print(f"  filter_mode = {filter_mode}, filter_threshold = {filter_threshold}", flush=True)
    print(f"  eim_span = {eim_span}", flush=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 读取 mzML
    rt_raw, rt_to_ions, ccs_spec = read_mzml_and_extract_ions(input_file, peak_intensity_thres, rt_range_max)
    np.save(os.path.join(output_dir, 'rt_to_ions.npy'), rt_to_ions)
    np.save(os.path.join(output_dir, 'ccs_spec.npy'), ccs_spec)   # 保存 CCS 数组

    # 2. 筛选 RT 范围
    rt, index = filter_rt_range(rt_raw, rt_range_min, rt_range_max)
    np.save(os.path.join(output_dir, 'rt.npy'), rt)
    np.save(os.path.join(output_dir, 'index.npy'), index)

    # 3. 提取 marker 强度
    marker_intensity = extract_marker_intensity(rt, rt_to_ions, marker, mzppm, ccsppm)
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
    df_cell = pd.DataFrame({'cell-rt': [float(x[0]) for x in cell_pulse],
                            'marker': [max(x) for x in cell_marker_vals]})
    df_cell.to_excel(os.path.join(output_dir, 'cell_marker.xlsx'), index=False)

    # 5. 特征峰聚类（TIMS 版本）
    peak_tar_filter = cluster_features_tims(ccs_spec, rt_to_ions, cell_pulse, mzppm, eim_span)
    np.save(os.path.join(output_dir, 'peak_tar_filter.npy'), peak_tar_filter)

    # 6. 构建 EIC 矩阵
    eic_matrix = extract_eic_matrix_tims(rt, rt_to_ions, peak_tar_filter, mzppm, ccsppm)
    np.save(os.path.join(output_dir, 'c_peak_eic_intensity.npy'), eic_matrix)

    # 7. 提取细胞脉冲内峰最大值
    c_peak = extract_cell_peak_max_tims(rt, cell_pulse, eic_matrix, peak_tar_filter)
    np.save(os.path.join(output_dir, 'c_peak.npy'), c_peak)

    # 8. 创建分组并提取子集
    path_p, cell_pulse_p, c_peak_p, eic_matrix_p, RT_sub = create_groups_and_subset(
        output_dir, rt, cell_pulse, c_peak, eic_matrix, bg_marker, bg_intensity_max, mzppm, cell_spec_lim, marker)

    # 9. 平均强度和 RSD
    aver_intensity, rsd = compute_average_intensity_and_rsd(c_peak_p)
    np.save(os.path.join(path_p, 'c_peak_aver_intensity.npy'), aver_intensity)
    np.save(os.path.join(path_p, 'c_peak_rsd.npy'), rsd)

    # 10. 背景统计（TIMS 版本）
    blank_aver, blank_std, blank_rsd, _ = compute_background_stats_tims(
        path_p, bg_marker, bg_intensity_max, mzppm, ccsppm,
        peak_tar_filter, eic_matrix_p, cell_pulse_p, RT_sub)
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
    fp_props, trash_props, sn, pulse_ratio, rsd_subs, total_score = compute_scores_and_filter_tims(
        peak_tar_filter, cell_occur, peak_pulse, aver_intensity,
        blank_aver, blank_std, blank_rsd, rsd, len(cell_pulse_p),
        filter_mode, filter_threshold)

    np.save(os.path.join(path_p, 'c_peak_sn.npy'), sn)
    np.save(os.path.join(path_p, 'c_peak_pulse_ratio.npy'), pulse_ratio)
    np.save(os.path.join(path_p, 'c_peak_rsd_subs.npy'), rsd_subs)

    # Excel 保存
    df_fp = pd.DataFrame(fp_props)
    df_fp.to_excel(os.path.join(path_p, 'c_peak_fp.xlsx'), index=False)
    df_trash = pd.DataFrame(trash_props)
    df_trash.to_excel(os.path.join(path_p, 'c_peak_trash.xlsx'), index=False)
    df_all = pd.DataFrame({
        'm/z': peak_tar_filter[:, 1],
        'ccs': peak_tar_filter[:, 0],
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
    np.save(os.path.join(path_p, 'c_peak_fp_ccs.npy'), np.column_stack((fp_props['ccs'], fp_props['mz'])))
    np.save(os.path.join(path_p, 'c_peak_trash_ccs.npy'), np.column_stack((trash_props['ccs'], trash_props['mz'])))

    # 14. 绘图
    plot_eics_tims(path_p, RT_sub, eic_matrix_p, peak_tar_filter, marker,
                   mzppm, ccsppm, fp_props, trash_props, plot_option,
                   plot_rt_min, plot_rt_max)
    # 最终汇总
    print("\n========== END ==========", flush=True)
    print(f"cell pulse: {len(cell_pulse)}", flush=True)
    print(f"peak list length: {len(peak_tar_filter)}", flush=True)
    print(f"cell-specific peak: {len(fp_props['mz'])}", flush=True)
    print(f"trash peak: {len(trash_props['mz'])}", flush=True)
    print(f"All results saved in: {output_dir}", flush=True)
    print(f"All subgroup results saved in: {path_p}", flush=True)


# In[ ]:


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

