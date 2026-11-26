#!/usr/bin/env python
# coding: utf-8
import sys
import os
import json
import math
import bisect
import warnings
import atexit

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pyteomics import mzml

warnings.filterwarnings('ignore')


# ------------------------------------------------------------
# 日志：写入 output_dir/output.txt，同时打印到 stdout
# ------------------------------------------------------------
def setup_logging(output_dir, also_print=True):
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    log_file = open(os.path.join(output_dir, 'output.txt'), 'w',
                    encoding='utf-8', buffering=1)

    class Tee:
        def __init__(self, *files):
            self.files = files

        def write(self, data):
            for f in self.files:
                try:
                    f.write(data)
                    f.flush()
                except (ValueError, OSError):
                    pass

        def flush(self):
            for f in self.files:
                try:
                    f.flush()
                except (ValueError, OSError):
                    pass

        def reconfigure(self, **kwargs):
            pass

    if also_print:
        sys.stdout = Tee(orig_stdout, log_file)
        sys.stderr = Tee(orig_stderr, log_file)
    else:
        sys.stdout = log_file
        sys.stderr = log_file

    def cleanup():
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        log_file.close()

    atexit.register(cleanup)


# ------------------------------------------------------------
# 绘图函数（no_tims 版本）
# ------------------------------------------------------------
def plot_eics(path_p, RT_sub, eic_matrix_p, c_peak_mz_all, marker, ppm,
              fp_mz, fp_sn, trash_mz, trash_score, plot_option,
              rt_min_plot, rt_max_plot):
    print(f"[plot] option={plot_option}, RT=[{rt_min_plot}, {rt_max_plot}]", flush=True)
    if plot_option == 'none':
        return

    # 找到 marker 对应的 EIC
    marker_idx = None
    for idx, mz_val in enumerate(c_peak_mz_all):
        if abs(mz_val - marker) < marker * ppm:
            marker_idx = idx
            break
    if marker_idx is None:
        print("     No marker EIC detected, skip", flush=True)
        return
    eic_marker = eic_matrix_p[marker_idx]

    def _plot_one(mz_val, sn_or_score, subdir, prefix):
        pos = np.where(np.abs(c_peak_mz_all - mz_val) < 1e-6)[0]
        if len(pos) == 0:
            return False
        intensity = eic_matrix_p[pos[0]]
        plt.figure(figsize=(10, 4))
        plt.subplot(2, 1, 1)
        plt.plot(RT_sub, eic_marker, label='marker')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Int.')
        plt.legend()
        plt.subplot(2, 1, 2)
        plt.plot(RT_sub, intensity, label=f'm/z {mz_val:.4f}')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Int.')
        plt.legend()
        plt.savefig(os.path.join(subdir, f'{prefix}{sn_or_score:.2f} mz{mz_val:.4f}.png'))
        plt.close()
        return True

    if plot_option in ('fp', 'all'):
        path_fp = os.path.join(path_p, 'plot_cell-specific')
        os.makedirs(path_fp, exist_ok=True)
        cnt = 0
        for idx, mz_val in enumerate(fp_mz):
            if _plot_one(mz_val, fp_sn[idx], path_fp, 'sn'):
                cnt += 1
        print(f"     Drawn {cnt} FP EIC plots", flush=True)

    if plot_option in ('trash', 'all'):
        path_trash = os.path.join(path_p, 'plot_trash')
        os.makedirs(path_trash, exist_ok=True)
        cnt = 0
        for idx, mz_val in enumerate(trash_mz):
            if _plot_one(mz_val, trash_score[idx], path_trash, 'score'):
                cnt += 1
        print(f"     Drawn {cnt} trash EIC plots", flush=True)


# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Error: no params file path given", file=sys.stderr, flush=True)
        sys.exit(1)

    with open(sys.argv[1], 'r', encoding='utf-8') as f:
        params = json.load(f)

    # ---------- 参数 ----------
    peak_intensity_thres = float(params['peak_intensity_thres'])
    ppm                  = float(params['ppm'])
    rt_range_min         = float(params['rt_range_min'])
    rt_range_max         = float(params['rt_range_max'])
    cell_eic_span        = int(params['cell_eic_span'])
    cell_eic_span_lim    = int(params['cell_eic_span_lim']) # cell pulse 最短长度
    cell_occur_thres     = float(params.get('cell_occurrence_threshold', 0.05))
    marker               = float(params['marker'])
    intensity_thres      = float(params['intensity_thres'])
    bg_marker            = np.array(params['bg_marker'], dtype=float)
    bg_intensity_max     = np.array(params['bg_intensity_max'], dtype=float)

    input_file           = params['input_file']
    output_dir           = params['output_dir']
    plot_option          = params.get('plot_option', 'none')
    plot_rt_min          = float(params.get('plot_rt_min', 5.0))
    plot_rt_max          = float(params.get('plot_rt_max', 10.0))
    filter_threshold     = float(params.get('filter_threshold', 10.0))
    spec_num_lim       = float(params.get('spec_num_lim', 1.0))

    # notebook 中的硬编码
    delta_spec     = 5      # 最大值范围
    

    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir, also_print=True)

    print("\n========== Starting MPCell (no_IM) ==========", flush=True)
    print(f"peak_intensity_thres = {peak_intensity_thres}", flush=True)
    print(f"ppm = {ppm}", flush=True)
    print(f"rt_range = [{rt_range_min}, {rt_range_max}]", flush=True)
    print(f"cell_eic_span = {cell_eic_span}, cell_eic_span_lim = {cell_eic_span_lim}", flush=True)
    print(f"cell_occur_thres = {cell_occur_thres}", flush=True)
    print(f"marker = {marker}, intensity_thres = {intensity_thres}", flush=True)
    print(f"bg_marker = {bg_marker}", flush=True)
    print(f"bg_intensity_max = {bg_intensity_max}", flush=True)
    print(f"input_file = {input_file}", flush=True)
    print(f"output_dir = {output_dir}", flush=True)
    print(f"plot_option = {plot_option}, plot_rt = [{plot_rt_min}, {plot_rt_max}]", flush=True)

    # --------------------------------------------------------
    # [1/14] 读取 mzML，提取 MS1 峰
    # --------------------------------------------------------
    print("[1/14] Read mzML file", flush=True)
    reader = mzml.read(input_file, use_index=False)
    rt = []
    rt_to_ions = {}
    rt_unit_detected = None
    for spectrum in reader:
        if spectrum.get('ms level') == 1:
            # 读取 scan start time 及其单位（mzML 规范单位为 minute 或 second）
            try:
                scan = spectrum['scanList']['scan'][0]
                rt_obj = scan['scan start time']
                rt_val = float(rt_obj)
                unit = getattr(rt_obj, 'unit_info', None)
            except (KeyError, IndexError, TypeError):
                continue

            # 记录首次出现的单位，供日志输出
            if rt_unit_detected is None and unit is not None:
                rt_unit_detected = unit
                print(f"     mzML RT unit detected: '{unit}'", flush=True)

            # 统一转换为秒
            if unit == 'minute':
                rt0 = rt_val * 60.0
            elif unit == 'second' or unit is None:
                rt0 = rt_val
            else:
                # 非预期单位（极少见），按秒处理并警告一次
                if rt_unit_detected != 'unknown_warned':
                    print(f"     Warning: unknown RT unit '{unit}', treated as second", flush=True)
                    rt_unit_detected = 'unknown_warned'
                rt0 = rt_val

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

    np.save(os.path.join(output_dir, 'rt_to_ions.npy'), rt_to_ions)

    rt = np.array(rt)
    index = np.array(range(len(rt)))
    mask_rt = (rt > rt_range_min) & (rt < rt_range_max)
    rt = rt[mask_rt]
    index = index[mask_rt]

    np.save(os.path.join(output_dir, 'rt.npy'), rt)
    np.save(os.path.join(output_dir, 'index.npy'), index)
    if rt_unit_detected is None:
        print("     Warning: no RT unit detected in mzML, treated as second", flush=True)
    print(f"     effective RT points: {len(rt)} (RT unit: second)", flush=True)
    print("[1/14] Finish", flush=True)

    if len(rt) == 0:
        print("Error: no RT points in specified range", file=sys.stderr, flush=True)
        sys.exit(1)

    # --------------------------------------------------------
    # [2/14] marker EIC
    # --------------------------------------------------------
    print(f"[2/14] Extract marker EIC (m/z={marker})", flush=True)
    marker_intensity = np.zeros(len(rt))
    for i, rt_val in enumerate(rt):
        if rt_val not in rt_to_ions:
            continue
        mz_arr, int_arr = rt_to_ions[rt_val]
        tol = marker * ppm
        m = abs(mz_arr - marker) < tol
        if np.any(m):
            marker_intensity[i] = int_arr[m][0]

    pd.DataFrame({'rt': rt, 'marker': marker_intensity}) \
        .to_excel(os.path.join(output_dir, 'marker_intensity.xlsx'), index=False)
    print("[2/14] Finish", flush=True)

    # --------------------------------------------------------
    # [3/14] cell pulse 检测
    # --------------------------------------------------------
    print(f"[3/14] Cell pulse detection (intensity_thres={intensity_thres})", flush=True)
    # spec_range 取决于 marker 的 m/z 大小
    if marker > 600:
        if cell_eic_span % 2 == 0:
            spec_range = list(range(int(-cell_eic_span / 2), int(cell_eic_span / 2)))
        else:
            spec_range = list(range(int(-cell_eic_span / 2), int(cell_eic_span / 2) + 1))
    else:
        if cell_eic_span % 2 == 0:
            spec_range = list(range(int(-cell_eic_span / 2) + 1, int(cell_eic_span / 2) + 1))
        else:
            spec_range = list(range(int(-cell_eic_span / 2), int(cell_eic_span / 2) + 1))

    mask = marker_intensity > intensity_thres
    mask_idx = np.where(mask)[0]
    if len(mask_idx) == 0:
        print("Warning: no marker signal above threshold, terminated", flush=True)
        sys.exit(1)

    cell_index = index[mask_idx]
    diff_idx = np.diff(cell_index)
    break_points = np.where(diff_idx != 1)[0] + 1
    groups = np.split(mask_idx, break_points)
    groups = [g for g in groups if len(g) >= cell_eic_span_lim]

    cell_pulse = []
    cell_marker = []
    for g in groups:
        group_intensity = marker_intensity[g]
        local_max_idx = np.argmax(group_intensity)
        global_max_idx = g[local_max_idx]
        start = global_max_idx + spec_range[0]
        end = global_max_idx + spec_range[-1]
        if start > 0 and end < len(rt):
            pulse_rts = rt[global_max_idx + np.array(spec_range)].tolist()
            cell_pulse.append(pulse_rts)
            pulse_marker = marker_intensity[global_max_idx + np.array(spec_range)].tolist()
            cell_marker.append(pulse_marker)

    # 相邻脉冲去重（notebook 中的 "叠脉冲" 过滤）
    cell_pulse_filter = []
    cell_marker_filter = []
    for idx in range(len(cell_pulse)):
        if idx == 0:
            # notebook 只处理 idx==1 和 idx>=2，因此保留 cell_pulse[0] 会被下面 idx==1 分支覆盖
            pass
        elif idx == 1:
            if np.median(cell_pulse[1]) >= cell_pulse[0][-1]:
                cell_pulse_filter.append(cell_pulse[0])
                cell_pulse_filter.append(cell_pulse[1])
                cell_marker_filter.append(cell_marker[0])
                cell_marker_filter.append(cell_marker[1])
            elif np.max(cell_marker[1]) > np.max(cell_marker[0]):
                cell_pulse_filter.append(cell_pulse[1])
                cell_marker_filter.append(cell_marker[1])
            else:
                cell_pulse_filter.append(cell_pulse[0])
                cell_marker_filter.append(cell_marker[0])
        else:  # idx >= 2
            if np.median(cell_pulse[idx]) >= cell_pulse_filter[-1][-1]:
                cell_pulse_filter.append(cell_pulse[idx])
                cell_marker_filter.append(cell_marker[idx])
            elif np.max(cell_marker[idx]) > np.max(cell_marker_filter[-1]):
                cell_pulse_filter.pop()
                cell_marker_filter.pop()
                cell_pulse_filter.append(cell_pulse[idx])
                cell_marker_filter.append(cell_marker[idx])
            elif np.max(cell_marker[idx]) == np.max(cell_marker_filter[-1]):
                max_idx = cell_marker[idx].index(np.max(cell_marker[idx]))
                max_idx0 = cell_marker_filter[-1].index(np.max(cell_marker_filter[-1]))
                if abs(len(cell_marker[idx]) - max_idx) < abs(len(cell_marker[idx]) - max_idx0):
                    cell_pulse_filter.pop()
                    cell_marker_filter.pop()
                    cell_pulse_filter.append(cell_pulse[idx])
                    cell_marker_filter.append(cell_marker[idx])

    cell_pulse = cell_pulse_filter
    cell_marker = cell_marker_filter

    if len(cell_pulse) == 0:
        print("Warning: no cell pulse after filtering, terminated", flush=True)
        sys.exit(1)

    np.save(os.path.join(output_dir, 'cell_pulse.npy'), np.array(cell_pulse, dtype=object))
    pd.DataFrame({
        'cell-rt': [float(x[0]) for x in cell_pulse],
        'marker': [max(x) for x in cell_marker]
    }).to_excel(os.path.join(output_dir, 'cell_marker.xlsx'), index=False)
    print(f"     cell pulses = {len(cell_pulse)}", flush=True)
    print("[3/14] Finish", flush=True)

    # --------------------------------------------------------
    # [4/14] 聚类提取 feature m/z
    # --------------------------------------------------------
    print(f"[4/14] Peak clustering (ppm={ppm})", flush=True)
    cell_pulse_trunc = cell_pulse[:100] if len(cell_pulse) >= 100 else cell_pulse
    rt_list = []
    for pulse in cell_pulse_trunc:
        for v in pulse:
            rt_list.append(float(v))
    rt_list = sorted(set(rt_list))

    mz_list = []
    for rt_val in rt_list:
        if rt_val not in rt_to_ions:
            continue
        mz_arr, int_arr = rt_to_ions[rt_val]
        for idx, mz in enumerate(mz_arr):
            mz_list.append((mz, int_arr[idx]))
    mz_list.sort(key=lambda x: x[1], reverse=True)
    mz_list = np.array(mz_list)
    print(f"     initial mz_list length = {len(mz_list)}", flush=True)

    cycle = 0
    while True:
        clusters = []
        ref_mzs = []
        for (mz, intensity) in mz_list:
            pos = bisect.bisect_left(ref_mzs, mz)
            found = False
            candidate_int = -1
            candidate_idx = -1
            for idx in (pos, pos - 1):
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
                clusters.insert(pos, {'ref_mz': mz, 'points': [(mz, intensity)]})
                ref_mzs.insert(pos, mz)

        percent = 0.1 if cycle == 0 else 0
        mz_list_final = []
        for cluster in clusters:
            points = np.array(cluster['points'])
            if len(points) > percent * len(cell_pulse_trunc):
                plus_sum = np.sum(points[:, 0] * points[:, 1])
                intensity_sum = np.sum(points[:, 1])
                mz_tar_v = plus_sum / intensity_sum
                mz_list_final.append((mz_tar_v, intensity_sum))
        mz_list_final.sort(key=lambda x: x[1], reverse=True)
        cycle += 1
        print(f"     iter {cycle}: n_cluster={len(clusters)}, n_mz={len(mz_list_final)}", flush=True)
        if len(mz_list) == len(mz_list_final):
            break
        mz_list = np.array(mz_list_final)

    mz_tar = np.array([x[0] for x in mz_list_final])
    np.save(os.path.join(output_dir, 'mz_tar.npy'), mz_tar)
    print(f"     final mz_tar length = {len(mz_tar)}", flush=True)
    print("[4/14] Finish", flush=True)

    # --------------------------------------------------------
    # [5/14] EIC 矩阵
    # --------------------------------------------------------
    print(f"[5/14] Build EIC matrix (n_peak={len(mz_tar)}, n_time={len(rt)})", flush=True)
    sorted_idx = np.argsort(mz_tar)
    mz_tar_sorted = mz_tar[sorted_idx]
    orig_order = np.argsort(sorted_idx)

    n_times = len(rt)
    n_peaks = len(mz_tar)
    time_intensity = np.zeros((n_peaks, n_times), dtype=np.float32)
    for t_idx, rt0 in enumerate(rt):
        if rt0 not in rt_to_ions:
            continue
        mz_arr, int_arr = rt_to_ions[rt0]
        for peak_idx in range(n_peaks):
            mz_target = mz_tar_sorted[peak_idx]
            tol = mz_target * ppm
            left = bisect.bisect_left(mz_arr, mz_target - tol)
            right = bisect.bisect_right(mz_arr, mz_target + tol)
            if left < right:
                time_intensity[peak_idx, t_idx] = int_arr[left:right][0]
    eic_matrix = time_intensity[orig_order, :]
    np.save(os.path.join(output_dir, 'c_peak_eic_intensity.npy'), eic_matrix)
    print("[5/14] Finish", flush=True)

    # --------------------------------------------------------
    # [6/14] c_peak：每个 pulse 内的最大值
    # --------------------------------------------------------
    print("[6/14] c_peak max within pulses", flush=True)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(rt)}
    c_peak_max = []
    for each in range(n_peaks):
        intensities = eic_matrix[each, :]
        rt_max_list = []
        intensity_max_list = []
        for pulse in cell_pulse:
            pulse_use = pulse
            if cell_eic_span > delta_spec * 2:
                mark = int(cell_eic_span / 2) - delta_spec
                if mark > 0:
                    pulse_use = pulse[mark:-mark]
                elif mark == 0:
                    pulse_use = pulse[1:-1]
            max_intensity = -1
            max_rt = None
            for rt0 in pulse_use:
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
                intensity_max_list.append(0.0)
        c_peak_max.append([mz_tar[each], rt_max_list, intensity_max_list])

    c_peak_arr = np.array(c_peak_max, dtype=object)
    np.save(os.path.join(output_dir, 'c_peak.npy'), c_peak_arr)
    print("[6/14] Finish", flush=True)

    # --------------------------------------------------------
    # [7/14] RT 分组
    # --------------------------------------------------------
    print("[7/14] RT grouping (RT_del)", flush=True)
    rt_list_all = rt.tolist()
    rt_del_val = int(rt_list_all[-1] // 4)
    if rt_del_val > 300:
        rt_del_val = 300
    RT_5min = [r for r in rt_list_all if r > rt_del_val]

    path_group = os.path.join(output_dir, 'grouping')
    os.makedirs(path_group, exist_ok=True)
    np.save(path_group+'/RT_del.npy',RT_5min)

    part = 'RT_del'
    path_p = os.path.join(path_group, part)
    os.makedirs(path_p, exist_ok=True)
    RT = np.load(os.path.join(path_group, part + '.npy')).tolist()
    np.save(os.path.join(path_p, 'RT.npy'), RT)
    print("[7/14] Finish", flush=True)

    # --------------------------------------------------------
    # [8/14] 分组取子集
    # --------------------------------------------------------
    print("[8/14] Subset extraction per group", flush=True)
    RT_arr = np.array(RT)
    rt_min, rt_max = RT_arr[0], RT_arr[-1]

    start_times = np.array([p[0] for p in cell_pulse])
    end_times = np.array([p[-1] for p in cell_pulse])
    keep = (start_times >= rt_min) & (end_times <= rt_max)
    cell_pulse_indices = np.where(keep)[0]
    if len(cell_pulse_indices) == 0:
        print("Warning: no cell pulse in this RT group, terminated", flush=True)
        sys.exit(1)

    cell_pulse_p = [cell_pulse[i] for i in cell_pulse_indices]
    start_cell_ind = cell_pulse_indices[0]
    end_cell_ind = cell_pulse_indices[-1]

    c_peak_p = []
    for i in range(len(c_peak_arr)):
        mz_v = c_peak_arr[i][0]
        rt_sub = c_peak_arr[i][1][start_cell_ind:end_cell_ind + 1]
        int_sub = c_peak_arr[i][2][start_cell_ind:end_cell_ind + 1]
        c_peak_p.append([mz_v, rt_sub, int_sub])

    rt_mask = (rt >= rt_min) & (rt <= rt_max)
    RT_indices = np.where(rt_mask)[0]
    start_rt_ind, end_rt_ind = RT_indices[0], RT_indices[-1]
    eic_matrix_p = eic_matrix[:, start_rt_ind:end_rt_ind + 1]

    np.save(os.path.join(path_p, 'cell_pulse.npy'), np.array(cell_pulse_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak.npy'), np.array(c_peak_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak_eic_intensity.npy'), eic_matrix_p)
    print(f"     pulses={len(cell_pulse_p)}, peaks={len(c_peak_p)}", flush=True)
    print("[8/14] Finish", flush=True)

    # --------------------------------------------------------
    # [9/14] c_peak pulse 内平均强度与 RSD
    # --------------------------------------------------------
    print("[9/14] c_peak averaged intensity / RSD", flush=True)
    n_peaks_p = len(c_peak_p)
    aver_intensity = np.zeros(n_peaks_p, dtype=float)
    rsd = np.zeros(n_peaks_p, dtype=float)
    c_peak_cell_dynamic_range = np.zeros(n_peaks_p)
    for i in range(n_peaks_p):
        intensity = c_peak_p[i][2]
        if not isinstance(intensity, np.ndarray):
            intensity = np.array(intensity)
        non_zero = intensity[intensity != 0]
        if len(non_zero) != 0:
            c_peak_cell_dynamic_range[i] = np.max(non_zero) - np.min(non_zero)
        sorted_vals = np.sort(non_zero)
        t = int(0.1 * len(sorted_vals))
        trimmed = sorted_vals[t:-t] if t > 0 else sorted_vals
        if len(trimmed) != 0:
            aver_intensity[i] = np.mean(trimmed)
        if aver_intensity[i] != 0:
            rsd[i] = np.std(trimmed) / aver_intensity[i]
    np.save(os.path.join(path_p, 'c_peak_aver_intensity.npy'), aver_intensity)
    np.save(os.path.join(path_p, 'c_peak_rsd.npy'), rsd)
    np.save(os.path.join(path_p, 'c_peak_cell_dynamic_range.npy'), c_peak_cell_dynamic_range)
    print("[9/14] Finish", flush=True)

    # --------------------------------------------------------
    # [10/14] 背景统计
    # --------------------------------------------------------
    print("[10/14] Background statistics", flush=True)
    cell_pulse_rt = []
    for pulse in cell_pulse_p:
        for rt0 in pulse:
            cell_pulse_rt.append(rt0)
    cell_pulse_rt = sorted(set(cell_pulse_rt))

    rt_to_index_p = {rt0: idx0 for idx0, rt0 in enumerate(RT)}
    cell_pulse_idx = set(rt_to_index_p[rt0] for rt0 in cell_pulse_rt if rt0 in rt_to_index_p)
    bg_nocell_idx = np.array([idx for idx in range(len(RT)) if idx not in cell_pulse_idx])

    n_marker = len(bg_marker)
    c_peak_blank_index = np.array([], dtype=int)
    for count in range(n_marker):
        bg_marker_each = bg_marker[count]
        bg_intensity_max_each = bg_intensity_max[count]
        m = np.abs(mz_tar - bg_marker_each) / bg_marker_each < ppm
        peak_indices = np.where(m)[0]
        if len(peak_indices) == 0:
            print(f"     Warning: no peak matches background marker {bg_marker_each}", flush=True)
            continue
        eic_intensity_marker = eic_matrix_p[peak_indices[0]]
        bg_intensity_idx = np.where(eic_intensity_marker <= bg_intensity_max_each)[0]
        if len(c_peak_blank_index) == 0:
            c_peak_blank_index = bg_intensity_idx
        else:
            c_peak_blank_index = np.intersect1d(c_peak_blank_index, bg_intensity_idx)

    c_peak_blank_index = np.intersect1d(bg_nocell_idx, c_peak_blank_index)
    print(f"     BG time points = {len(c_peak_blank_index)} / {len(RT)}", flush=True)

    spectra_num = len(eic_matrix_p[0])
    if spectra_num > 0 and len(c_peak_blank_index) / spectra_num < 0.33:
        print("     Warning: background spectra < 1/3 total spectra!", flush=True)
    if len(c_peak_blank_index) == 0:
        print("Error: no background spectra defined", file=sys.stderr, flush=True)
        sys.exit(1)

    n_peaks_2 = eic_matrix_p.shape[0]
    c_peak_blank_aver = np.zeros(n_peaks_2)
    c_peak_blank_rsd = np.zeros(n_peaks_2)
    c_peak_blank_std = np.zeros(n_peaks_2)
    c_peak_blank_dynamic_range = np.zeros(n_peaks_2)
    c_peak_blank_aver_raw = np.zeros(n_peaks_2)
    c_peak_blank_rsd_raw = np.zeros(n_peaks_2)
    c_peak_blank_std_raw = np.zeros(n_peaks_2)

    for i in range(n_peaks_2):
        blank_intensity = eic_matrix_p[i, c_peak_blank_index]
        if len(blank_intensity) > 0:
            t = int(0.1 * len(blank_intensity))
            trimmed_raw = blank_intensity[t:-t] if t > 0 else blank_intensity
            c_peak_blank_aver_raw[i] = np.mean(trimmed_raw)
            c_peak_blank_std_raw[i] = np.std(trimmed_raw)
            if c_peak_blank_aver_raw[i] != 0:
                c_peak_blank_rsd_raw[i] = c_peak_blank_std_raw[i] / c_peak_blank_aver_raw[i]

        non_zero = blank_intensity[blank_intensity != 0]
        if len(non_zero) != 0:
            sorted_vals = np.sort(non_zero)
            topn = math.floor(0.1 * len(sorted_vals))
            trimmed = sorted_vals[topn:-topn] if topn > 0 else sorted_vals
            c_peak_blank_aver[i] = np.mean(trimmed)
            c_peak_blank_std[i] = np.std(trimmed)
            if c_peak_blank_aver[i] != 0:
                c_peak_blank_rsd[i] = c_peak_blank_std[i] / c_peak_blank_aver[i]
            c_peak_blank_dynamic_range[i] = np.max(non_zero) - np.min(non_zero)

    np.save(os.path.join(path_p, 'c_peak_blank_aver_intensity.npy'), c_peak_blank_aver)
    np.save(os.path.join(path_p, 'c_peak_blank_rsd.npy'), c_peak_blank_rsd)
    np.save(os.path.join(path_p, 'c_peak_blank_std.npy'), c_peak_blank_std)
    np.save(os.path.join(path_p, 'c_peak_blank_dynamic_range.npy'), c_peak_blank_dynamic_range)
    np.save(os.path.join(path_p, 'c_peak_blank_aver_intensity_raw.npy'), c_peak_blank_aver_raw)
    np.save(os.path.join(path_p, 'c_peak_blank_rsd_raw.npy'), c_peak_blank_rsd_raw)
    np.save(os.path.join(path_p, 'c_peak_blank_std_raw.npy'), c_peak_blank_std_raw)
    print("[10/14] Finish", flush=True)

    # --------------------------------------------------------
    # [11/14] cell_occur
    # --------------------------------------------------------
    print("[11/14] cell_occur", flush=True)
    cell_occur = np.zeros(n_peaks_2, dtype=np.int64)
    cell_occur_raw = np.zeros(n_peaks_2, dtype=np.int64)

    for each in range(n_peaks_2):
        intensities = eic_matrix_p[each, :]
        blank_val = c_peak_blank_aver_raw[each]
        blank_std = c_peak_blank_std_raw[each]
        count = 0
        for pulse in cell_pulse_p:
            intensity_list = []
            for rt0 in pulse:
                idx = rt_to_index_p.get(rt0)
                if idx is not None:
                    intensity_list.append(intensities[idx])
            if len(intensity_list) > 0:
                idx1 = np.where(np.array(intensity_list) > blank_val + 3 * blank_std)[0]
                diff_idx = np.diff(idx1)
                if len(idx1) > 0 and len(np.where(diff_idx == 1)[0]) >= spec_num_lim - 1:
                    count += 1
        cell_occur_raw[each] = count

    for each in range(n_peaks_2):
        intensities = eic_matrix_p[each, :]
        blank_val = c_peak_blank_aver[each]
        blank_std = c_peak_blank_std[each]
        count = 0
        for pulse in cell_pulse_p:
            pulse_use = pulse
            if cell_eic_span > delta_spec * 2:
                mark = int(cell_eic_span / 2) - delta_spec
                if mark > 0:
                    pulse_use = pulse[mark:-mark]
                elif mark == 0:
                    pulse_use = pulse[1:-1]
            intensity_list = []
            for rt0 in pulse_use:
                idx = rt_to_index_p.get(rt0)
                if idx is not None:
                    intensity_list.append(intensities[idx])
            if len(intensity_list) > 0:
                idx1 = np.where(np.array(intensity_list) > blank_val + 3 * blank_std)[0]
                diff_idx = np.diff(idx1)
                if len(idx1) > 0 and len(np.where(diff_idx == 1)[0]) >= spec_num_lim - 1:
                    count += 1
        cell_occur[each] = count

    n_pulses = len(cell_pulse_p)
    cell_occur_ratio = cell_occur / n_pulses
    cell_occur_ratio_raw = cell_occur_raw / n_pulses
    np.save(os.path.join(path_p, 'cell_occur.npy'), cell_occur)
    np.save(os.path.join(path_p, 'cell_occur_ratio.npy'), cell_occur_ratio)
    np.save(os.path.join(path_p, 'cell_occur_ratio_raw.npy'), cell_occur_ratio_raw)
    print(f"     cell_occur mean = {np.mean(cell_occur):.2f}", flush=True)
    print("[11/14] Finish", flush=True)

    # --------------------------------------------------------
    # [12/14] c_peak_pulse
    # --------------------------------------------------------
    print("[12/14] c_peak_pulse", flush=True)
    c_peak_pulse1 = np.zeros(n_peaks_2, dtype=np.int64)
    c_peak_pulse2 = np.zeros(n_peaks_2, dtype=np.int64)
    for i in range(n_peaks_2):
        intensity = eic_matrix_p[i]
        threshold = c_peak_blank_aver_raw[i] + 3 * c_peak_blank_std_raw[i]
        pos = np.where(intensity > threshold)[0]
        if len(pos) > 0:
            diff = np.diff(pos)
            split_indices = np.where(diff != 1)[0] + 1
            groups_ = np.split(pos, split_indices)
            c_peak_pulse1[i] = sum(1 for g in groups_ if len(g) >= spec_num_lim)

    for i in range(n_peaks_2):
        intensity = eic_matrix_p[i]
        pos = np.where(intensity > 0)[0]
        if len(pos) > 0:
            diff = np.diff(pos)
            split_indices = np.where(diff != 1)[0] + 1
            groups_ = np.split(pos, split_indices)
            c_peak_pulse2[i] = sum(1 for g in groups_ if len(g) >= spec_num_lim)

    c_peak_pulse = np.maximum(c_peak_pulse1, c_peak_pulse2)
    np.save(os.path.join(path_p, 'c_peak_pulse.npy'), c_peak_pulse)
    print("[12/14] Finish", flush=True)

    # --------------------------------------------------------
    # [13/14] SN / pulse_ratio / rsd_subs
    # --------------------------------------------------------
    print("[13/14] SN, pulse ratio, rsd_subs", flush=True)
    pulse_ratio = np.zeros_like(c_peak_pulse, dtype=float)
    mask = c_peak_pulse != 0
    np.divide(cell_occur_raw, c_peak_pulse, out=pulse_ratio, where=mask)

    rsd_subs = rsd - c_peak_blank_rsd

    sn = np.zeros_like(aver_intensity, dtype=float)
    mask_rsd_non0 = c_peak_blank_rsd != 0
    sn[mask_rsd_non0] = (aver_intensity[mask_rsd_non0] - c_peak_blank_aver[mask_rsd_non0]) \
                        / c_peak_blank_std[mask_rsd_non0]
    mask_rsd0_blank_non0 = (~mask_rsd_non0) & (c_peak_blank_aver != 0)
    sn[mask_rsd0_blank_non0] = (aver_intensity[mask_rsd0_blank_non0]
                                - c_peak_blank_aver[mask_rsd0_blank_non0]) \
                               / c_peak_blank_aver[mask_rsd0_blank_non0]

    np.save(os.path.join(path_p, 'c_peak_sn.npy'), sn)
    np.save(os.path.join(path_p, 'c_peak_pulse_ratio.npy'), pulse_ratio)
    np.save(os.path.join(path_p, 'c_peak_rsd_subs.npy'), rsd_subs)
    print("[13/14] Finish", flush=True)

    # --------------------------------------------------------
    # [14/14] 多重滤除打分 + 筛选 + 绘图
    # --------------------------------------------------------
    print("[14/14] Multi-parameter scoring and filtering", flush=True)
    eff_cell = np.zeros(n_peaks_p)
    for idx, each in enumerate(c_peak_p):
        cell_intensity = np.array(each[2], dtype=float)
        if len(cell_intensity) > 0:
            eff_cell[idx] = len(cell_intensity[cell_intensity > 0]) / len(cell_intensity)

    # ---- 读取用户配置的评分参数（默认与原逻辑一致）----
    use_snr   = bool(params.get('use_snr',   True))
    use_match = bool(params.get('use_match', True))
    use_rsd   = bool(params.get('use_rsd',   True))

    snr_thres   = float(params.get('snr_thres',   3.5))
    snr_lower   = float(params.get('snr_lower',   2.0))
    match_thres = float(params.get('match_thres', 0.67))
    match_lower = float(params.get('match_lower', 0.3))
    rsd_thres   = float(params.get('rsd_thres',   0.2))
    rsd_lower   = float(params.get('rsd_lower',   0.1))

    print(f"     use_snr={use_snr} thres={snr_thres} lower={snr_lower}", flush=True)
    print(f"     use_match={use_match} thres={match_thres} lower={match_lower}", flush=True)
    print(f"     use_rsd={use_rsd} thres={rsd_thres} lower={rsd_lower}", flush=True)

    score = np.zeros(n_peaks_p)
    score_sn = np.zeros_like(sn)
    score_rsd = np.zeros_like(rsd_subs)
    score_match = np.zeros_like(pulse_ratio)

    # ---- SNR ----
    if use_snr:
        for idx in np.where((sn >= snr_thres) & (cell_occur_ratio > cell_occur_thres))[0]:
            score[idx] = 10
        for idx in np.where((sn < snr_thres) & (sn >= snr_lower)
                            & (cell_occur_ratio > cell_occur_thres))[0]:
            score_sn[idx] = 10 * (sn[idx] - snr_lower) / (snr_thres - snr_lower)

    # ---- RSD subs ----
    if use_rsd:
        for idx in np.where((rsd_subs >= rsd_thres) & (cell_occur_ratio > cell_occur_thres))[0]:
            score[idx] = 10
        for idx in np.where((rsd_subs < rsd_thres) & (rsd_subs >= rsd_lower)
                            & (cell_occur_ratio > cell_occur_thres))[0]:
            score_rsd[idx] = 10 * (rsd_subs[idx] - rsd_lower) / (rsd_thres - rsd_lower)

    # ---- Match ratio ----
    if use_match:
        for idx in np.where((pulse_ratio >= match_thres) & (eff_cell > 0.05))[0]:
            score[idx] = 10
        for idx in np.where((pulse_ratio < match_thres) & (pulse_ratio >= match_lower)
                            & (eff_cell > 0.05))[0]:
            score_match[idx] = 10 * (pulse_ratio[idx] - match_lower) / (match_thres - match_lower)

    for idx in range(len(score)):
        if score[idx] != 10:
            three = sorted([score_sn[idx], score_rsd[idx], score_match[idx]])
            score[idx] = sum(three)
            
    fp_mask = (score >= filter_threshold) & (cell_occur_ratio >= cell_occur_thres)    
    fp_mz = mz_tar[fp_mask]
    trash_mz = mz_tar[~fp_mask]
    fp_sn = sn[fp_mask]
    trash_score = score[~fp_mask]
    print(f"     Positive peaks = {len(fp_mz)}, Trash peaks = {len(trash_mz)}", flush=True)

    np.save(os.path.join(output_dir, 'cell_specific_peak.npy'), fp_mz)
    np.save(os.path.join(output_dir, 'cell_nonspecific_peak.npy'), trash_mz)
    np.save(os.path.join(path_p, 'cell_specific_peak.npy'), fp_mz)
    np.save(os.path.join(path_p, 'cell_nonspecific_peak.npy'), trash_mz)

    pd.DataFrame({
        'mz': mz_tar[fp_mask],
        'snr': sn[fp_mask],
        'match_ratio': pulse_ratio[fp_mask],
        'rsd_subs': rsd_subs[fp_mask],
        'cell_occur_ratio': cell_occur_ratio[fp_mask],
        'score': score[fp_mask],
    }).to_excel(os.path.join(path_p, 'cell_specific_peak.xlsx'), index=False)

    pd.DataFrame({
        'mz': mz_tar[~fp_mask],
        'snr': sn[~fp_mask],
        'match_ratio': pulse_ratio[~fp_mask],
        'rsd_subs': rsd_subs[~fp_mask],
        'cell_occur_ratio': cell_occur_ratio[~fp_mask],
        'score': score[~fp_mask],
    }).to_excel(os.path.join(path_p, 'cell_nonspecific_peak.xlsx'), index=False)

    pd.DataFrame({
        'm/z': mz_tar,
        'snr': sn,
        'match_ratio': pulse_ratio,
        'rsd_subs': rsd_subs,
        'cell_occur_ratio': cell_occur_ratio,
        'score': score,
    }).to_excel(os.path.join(path_p, 'all_candidate_peak.xlsx'), index=False)

    # 绘图
    plot_eics(path_p, np.array(RT), eic_matrix_p, mz_tar, marker, ppm,
              fp_mz, fp_sn, trash_mz, trash_score,
              plot_option, plot_rt_min, plot_rt_max)
    print("[14/14] Finish", flush=True)

    print("\n========== END ==========", flush=True)
    print(f"cell pulses      : {len(cell_pulse)}", flush=True)
    print(f"peak candidate   : {len(mz_tar)}", flush=True)
    print(f"cell-specific    : {len(fp_mz)}", flush=True)
    print(f"non-specific     : {len(trash_mz)}", flush=True)
    print(f"results saved to : {output_dir}", flush=True)
    print(f"subgroup saved to: {path_p}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)