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
# 绘图函数（tims 版本）
# ------------------------------------------------------------
def plot_eics_tims(path_p, RT_sub, eic_matrix_p, c_peak_mz_all, marker,
                   ppm_frac, ccs_range, fp_mz, fp_sn, trash_mz, trash_score,
                   plot_option, rt_min_plot, rt_max_plot):
    """
    c_peak_mz_all: object array，每个元素为 (mz, ccs)
    marker: (marker_mz, marker_ccs)
    ppm_frac: ppm 的小数形式（如 20ppm → 20e-6），# [FIX] 传入时已转好
    """
    print(f"[plot] option={plot_option}, RT=[{rt_min_plot}, {rt_max_plot}]", flush=True)
    if plot_option == 'none':
        return

    # 找到 marker 对应的 EIC
    marker_idx = None
    for idx, mz_ccs in enumerate(c_peak_mz_all):
        mz_val, ccs_val = mz_ccs
        # [FIX] 使用 ppm_frac 而不是 ppm * marker[0]
        if abs(mz_val - marker[0]) / marker[0] < ppm_frac and abs(ccs_val - marker[1]) < ccs_range:
            marker_idx = idx
            break
    if marker_idx is None:
        print("     No marker EIC detected, skip", flush=True)
        return
    eic_marker = eic_matrix_p[marker_idx]

    def _find_index(mz_ccs):
        mz_val, ccs_val = mz_ccs
        for i, p in enumerate(c_peak_mz_all):
            if abs(p[0] - mz_val)/mz_val < ppm_frac and abs(p[1] - ccs_val) < ccs_range:
                return i
        return None

    def _plot_one(mz_ccs, sn_or_score, subdir, prefix):
        idx = _find_index(mz_ccs)
        if idx is None:
            return False
        intensity = eic_matrix_p[idx]
        mz_val, ccs_val = mz_ccs
        plt.figure(figsize=(10, 4))
        plt.subplot(2, 1, 1)
        plt.plot(RT_sub, eic_marker, label='marker')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Int.')
        plt.legend()
        plt.subplot(2, 1, 2)
        plt.plot(RT_sub, intensity, label=f'm/z {mz_val:.4f} ccs {ccs_val:.3f}')
        plt.xlim(rt_min_plot, rt_max_plot)
        plt.xlabel('RT/s')
        plt.ylabel('Int.')
        plt.legend()
        plt.savefig(os.path.join(subdir, f'{prefix}{sn_or_score:.2f} mz{mz_val:.4f}_ccs{ccs_val:.3f}.png'))
        plt.close()
        return True

    if plot_option in ('fp', 'all'):
        path_fp = os.path.join(path_p, 'plot_cell-specific')
        os.makedirs(path_fp, exist_ok=True)
        cnt = 0
        for idx, mz_ccs in enumerate(fp_mz):
            if _plot_one(mz_ccs, fp_sn[idx], path_fp, 'score'):
                cnt += 1
        print(f"     Drawn {cnt} FP EIC plots", flush=True)

    if plot_option in ('trash', 'all'):
        path_trash = os.path.join(path_p, 'plot_trash')
        os.makedirs(path_trash, exist_ok=True)
        cnt = 0
        for idx, mz_ccs in enumerate(trash_mz):
            if _plot_one(mz_ccs, trash_score[idx], path_trash, 'score'):
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

    # ---------- 参数（全部 .get 容错）----------
    input_dir            = params.get('d_folder', '')
    output_dir           = params.get('output_dir', '')

    ppm                  = float(params.get('mz_tol_ppm', 20.0))
    ppm_frac             = ppm * 1e-6  # [FIX] ppm → 小数

    mobility_range       = float(params.get('mobility_range', 0.05))
    ccs_range            = float(params.get('ccs_range', mobility_range))

    rt_range_min         = float(params.get('rt_range_min', 0.0))
    rt_range_max         = float(params.get('rt_range_max', 1e9))

    cell_eic_span        = int(params.get('cell_eic_span', 5))
    cell_eic_span_lim    = int(params.get('cell_eic_span_lim', 1))
    spec_num_lim         = int(params.get('spec_num_lim', 1))
    cell_occur_thres     = float(params.get('cell_occurrence_threshold', 0.05))

    marker_mz            = float(params.get('marker_mz', 760.5845))
    marker_ccs           = float(params.get('marker_ccs',
                                            params.get('marker_mobility', 1.402)))
    intensity_thres      = float(params.get('int_lower', 60000.0))

    bg_marker_raw        = params.get('bg_marker', [])
    bg_intensity_max_raw = params.get('bg_intensity_max', [])
    if len(bg_marker_raw) > 0:
        bg_marker = np.array(bg_marker_raw, dtype=float).reshape(-1, 2)
    else:
        bg_marker = np.zeros((0, 2), dtype=float)
    bg_intensity_max = np.array(bg_intensity_max_raw, dtype=float) \
                       if len(bg_intensity_max_raw) > 0 else np.zeros((0,), dtype=float)

    plot_option          = params.get('plot_option', 'none')
    plot_rt_min          = float(params.get('plot_rt_min', 5.0))*60
    plot_rt_max          = float(params.get('plot_rt_max', 10.0))*60
    filter_threshold     = float(params.get('filter_threshold', 10.0))
    use_snr   = bool(params.get('use_snr',   True))
    use_match = bool(params.get('use_match', True))
    use_rsd   = bool(params.get('use_rsd',   True))

    snr_thres   = float(params.get('snr_thres',   3.5))
    snr_lower   = float(params.get('snr_lower',   2.0))
    match_thres = float(params.get('match_thres', 0.67))
    match_lower = float(params.get('match_lower', 0.3))
    rsd_thres   = float(params.get('rsd_thres',   0.2))
    rsd_lower   = float(params.get('rsd_lower',   0.1))

    delta_spec = 5  # 与 notebook 一致

    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir, also_print=True)

    # [FIX] d_folder 容错：如果根目录没有 npy，尝试 results 子目录
    if input_dir and not os.path.exists(
            os.path.join(input_dir, 'all_features_EIC_full_intensity.npy')):
        alt = os.path.join(input_dir, 'results')
        if os.path.exists(
                os.path.join(alt, 'all_features_EIC_full_intensity.npy')):
            print(f"[FIX] d_folder 自动切换到: {alt}", flush=True)
            input_dir = alt

    print("\n========== Starting MPCell (tims) ==========", flush=True)
    print(f"input_dir = {input_dir}", flush=True)
    print(f"output_dir = {output_dir}", flush=True)
    print(f"ppm = {ppm} (fraction={ppm_frac})", flush=True)
    print(f"rt_range = [{rt_range_min}, {rt_range_max}]", flush=True)
    print(f"cell_eic_span = {cell_eic_span}, cell_eic_span_lim = {cell_eic_span_lim}, spec_num_lim = {spec_num_lim}", flush=True)
    print(f"cell_occur_thres = {cell_occur_thres}", flush=True)
    print(f"marker = ({marker_mz}, {marker_ccs}), ccs_range = {ccs_range}", flush=True)
    print(f"intensity_thres = {intensity_thres}", flush=True)
    print(f"bg_marker = {bg_marker}", flush=True)
    print(f"bg_intensity_max = {bg_intensity_max}", flush=True)
    print(f"plot_option = {plot_option}, plot_rt = [{plot_rt_min}, {plot_rt_max}]", flush=True)
    print(f"use_snr={use_snr}  thres={snr_thres}  lower={snr_lower}", flush=True)
    print(f"use_match={use_match} thres={match_thres} lower={match_lower}", flush=True)
    print(f"use_rsd={use_rsd}  thres={rsd_thres}  lower={rsd_lower}", flush=True)

    # --------------------------------------------------------
    # [1/14] 读取已提取的 EIC 全矩阵和 feature 表
    # --------------------------------------------------------
    print("[1/14] Load EIC matrix and feature table", flush=True)
    eic_path = os.path.join(input_dir, 'all_features_EIC_full_intensity.npy')
    rt_path = os.path.join(input_dir, 'all_features_EIC_full_retention_time.npy')
    feat_path = os.path.join(input_dir, '01_feature_table.csv')
    for p in (eic_path, rt_path, feat_path):
        if not os.path.exists(p):
            print(f"Error: required file not found: {p}", file=sys.stderr, flush=True)
            sys.exit(1)

    eic = np.load(eic_path)
    rt = np.load(rt_path)
    df = pd.read_csv(feat_path)

    index = np.array(range(len(rt)))
    mask_rt = (rt > rt_range_min) & (rt < rt_range_max)
    rt = rt[mask_rt]
    index = index[mask_rt]

    c_peak_eic = []
    for peak in eic:
        c_peak_eic.append(peak[mask_rt])
    c_peak_eic = np.array(c_peak_eic, dtype=np.float32)

    mz_arr = df['mz'].values
    mob_arr = df['mobility'].values
    mz_tar = [(float(mz_arr[i]), float(mob_arr[i])) for i in range(len(mz_arr))]

    np.save(os.path.join(output_dir, 'rt.npy'), rt)
    np.save(os.path.join(output_dir, 'index.npy'), index)
    np.save(os.path.join(output_dir, 'c_peak_eic_intensity.npy'), c_peak_eic)
    np.save(os.path.join(output_dir, 'mz_tar.npy'), np.array(mz_tar, dtype=object))
    print(f"     effective RT points: {len(rt)}", flush=True)
    print(f"     number of peaks: {len(mz_tar)}", flush=True)
    print("[1/14] Finish", flush=True)

    if len(rt) == 0:
        print("Error: no RT points in specified range", file=sys.stderr, flush=True)
        sys.exit(1)

    # --------------------------------------------------------
    # [2/14] marker EIC
    # --------------------------------------------------------
    print(f"[2/14] Extract marker EIC (m/z={marker_mz}, ccs={marker_ccs})", flush=True)
    marker_intensity = None
    for idx, (mz_val, ccs_val) in enumerate(mz_tar):
        # [FIX] 使用 ppm_frac
        if abs(mz_val - marker_mz) / marker_mz < ppm_frac and abs(ccs_val - marker_ccs) < ccs_range:
            marker_intensity = c_peak_eic[idx]
            break
    if marker_intensity is None:
        print("Error: marker not found in feature table", file=sys.stderr, flush=True)
        sys.exit(1)

    pd.DataFrame({'rt': rt, 'marker': marker_intensity}) \
        .to_excel(os.path.join(output_dir, 'marker_intensity.xlsx'), index=False)
    print("[2/14] Finish", flush=True)

    # --------------------------------------------------------
    # [3/14] cell pulse 检测
    # --------------------------------------------------------
    print(f"[3/14] Cell pulse detection (intensity_thres={intensity_thres})", flush=True)
    if marker_mz > 600:
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

    cell_pulse_filter = []
    cell_marker_filter = []
    for idx in range(len(cell_pulse)):
        if idx == 0:
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
        else:
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
    # [4/14] c_peak
    # --------------------------------------------------------
    print("[4/14] c_peak max within pulses", flush=True)
    rt_to_index = {rt0: idx0 for idx0, rt0 in enumerate(rt)}
    c_peak_max = []
    for peak_idx in range(len(mz_tar)):
        intensities = c_peak_eic[peak_idx, :]
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
        c_peak_max.append([mz_tar[peak_idx], rt_max_list, intensity_max_list])

    c_peak_arr = np.array(c_peak_max, dtype=object)
    np.save(os.path.join(output_dir, 'c_peak.npy'), c_peak_arr)
    print("[4/14] Finish", flush=True)

    # --------------------------------------------------------
    # [5/14] RT 分组
    # --------------------------------------------------------
    print("[5/14] RT grouping (RT_del)", flush=True)
    rt_list_all = rt.tolist()
    rt_del_val = int(rt_list_all[-1] // 4)
    if rt_del_val > 300:
        rt_del_val = 300
    RT_5min = [r for r in rt_list_all if r > rt_del_val]

    path_group = os.path.join(output_dir, 'grouping')
    os.makedirs(path_group, exist_ok=True)
    np.save(os.path.join(path_group, 'RT_del.npy'), RT_5min)

    part = 'RT_del'
    path_p = os.path.join(path_group, part)
    os.makedirs(path_p, exist_ok=True)
    RT = np.load(os.path.join(path_group, part + '.npy')).tolist()
    np.save(os.path.join(path_p, 'RT.npy'), RT)
    print("[5/14] Finish", flush=True)

    # --------------------------------------------------------
    # [6/14] 分组取子集
    # --------------------------------------------------------
    print("[6/14] Subset extraction per group", flush=True)
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
        mz_ccs = c_peak_arr[i][0]
        rt_sub = c_peak_arr[i][1][start_cell_ind:end_cell_ind + 1]
        int_sub = c_peak_arr[i][2][start_cell_ind:end_cell_ind + 1]
        c_peak_p.append([mz_ccs, rt_sub, int_sub])

    rt_mask = (rt >= rt_min) & (rt <= rt_max)
    RT_indices = np.where(rt_mask)[0]
    start_rt_ind, end_rt_ind = RT_indices[0], RT_indices[-1]
    eic_matrix_p = c_peak_eic[:, start_rt_ind:end_rt_ind + 1]

    np.save(os.path.join(path_p, 'cell_pulse.npy'), np.array(cell_pulse_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak.npy'), np.array(c_peak_p, dtype=object))
    np.save(os.path.join(path_p, 'c_peak_eic_intensity.npy'), eic_matrix_p)
    print(f"     pulses={len(cell_pulse_p)}, peaks={len(c_peak_p)}", flush=True)
    print("[6/14] Finish", flush=True)

    # --------------------------------------------------------
    # [7/14] c_peak averaged intensity / RSD
    # --------------------------------------------------------
    print("[7/14] c_peak averaged intensity / RSD", flush=True)
    n_peaks_p = len(c_peak_p)
    aver_intensity = np.zeros(n_peaks_p, dtype=float)
    rsd = np.zeros(n_peaks_p, dtype=float)
    c_peak_cell_dynamic_range = np.zeros(n_peaks_p)
    for i in range(n_peaks_p):
        intensity = c_peak_p[i][2]
        if not isinstance(intensity, np.ndarray):
            intensity = np.array(intensity, dtype=float)
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
    print("[7/14] Finish", flush=True)

    # --------------------------------------------------------
    # [8/14] 背景统计
    # --------------------------------------------------------
    print("[8/14] Background statistics", flush=True)
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
    if n_marker > 0:
        for count in range(n_marker):
            bg_marker_each = bg_marker[count]
            bg_intensity_max_each = bg_intensity_max[count]
            found = False
            for idx, (mz_val, ccs_val) in enumerate(mz_tar):
                # [FIX] 使用 ppm_frac
                if abs(mz_val - bg_marker_each[0]) / bg_marker_each[0] < ppm_frac and \
                   abs(ccs_val - bg_marker_each[1]) < ccs_range:
                    eic_intensity_marker = eic_matrix_p[idx]
                    found = True
                    break
            if not found:
                print(f"     Warning: no peak matches background marker {bg_marker_each}", flush=True)
                continue
            bg_intensity_idx = np.where(eic_intensity_marker <= bg_intensity_max_each)[0]
            if len(c_peak_blank_index) == 0:
                c_peak_blank_index = bg_intensity_idx
            else:
                c_peak_blank_index = np.intersect1d(c_peak_blank_index, bg_intensity_idx)

        c_peak_blank_index = np.intersect1d(bg_nocell_idx, c_peak_blank_index)
    else:
        # [FIX] 无背景 marker 时：用所有非 cell pulse 点作为背景
        print("     No background marker given, using all non-cell time points", flush=True)
        c_peak_blank_index = bg_nocell_idx

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
    print("[8/14] Finish", flush=True)

    # --------------------------------------------------------
    # [9/14] cell_occur
    # --------------------------------------------------------
    print("[9/14] cell_occur", flush=True)
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
                arr = np.array(intensity_list)
                idx1 = np.where(arr > blank_val + 3 * blank_std)[0]
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
                arr = np.array(intensity_list)
                idx1 = np.where(arr > blank_val + 3 * blank_std)[0]
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
    print("[9/14] Finish", flush=True)

    # --------------------------------------------------------
    # [10/14] c_peak_pulse
    # --------------------------------------------------------
    print("[10/14] c_peak_pulse", flush=True)
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
    print("[10/14] Finish", flush=True)

    # --------------------------------------------------------
    # [11/14] SN / pulse_ratio / rsd_subs
    # --------------------------------------------------------
    print("[11/14] SN, pulse ratio, rsd_subs", flush=True)
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
    print("[11/14] Finish", flush=True)

    # --------------------------------------------------------
    # [12/14] 多重滤除打分 + 筛选 + 绘图
    # --------------------------------------------------------
    print("[12/14] Multi-parameter scoring and filtering", flush=True)
    eff_cell = np.zeros(n_peaks_p)
    for idx, each in enumerate(c_peak_p):
        cell_intensity = np.array(each[2], dtype=float)
        if len(cell_intensity) > 0:
            eff_cell[idx] = len(cell_intensity[cell_intensity > 0]) / len(cell_intensity)

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
    fp_mz = np.array([mz_tar[i] for i in np.where(fp_mask)[0]], dtype=object)
    trash_mz = np.array([mz_tar[i] for i in np.where(~fp_mask)[0]], dtype=object)
    fp_sn = sn[fp_mask]
    trash_score = score[~fp_mask]
    print(f"     Positive peaks = {len(fp_mz)}, Trash peaks = {len(trash_mz)}", flush=True)

    np.save(os.path.join(output_dir, 'cell_specific_peak.npy'), fp_mz)
    np.save(os.path.join(output_dir, 'cell_nonspecific_peak.npy'), trash_mz)
    np.save(os.path.join(path_p, 'cell_specific_peak.npy'), fp_mz)
    np.save(os.path.join(path_p, 'cell_nonspecific_peak.npy'), trash_mz)

    fp_mz_arr = np.array([list(x) for x in fp_mz]) if len(fp_mz) > 0 else np.empty((0, 2))
    trash_mz_arr = np.array([list(x) for x in trash_mz]) if len(trash_mz) > 0 else np.empty((0, 2))

    pd.DataFrame({
        'mz': fp_mz_arr[:, 0] if len(fp_mz_arr) > 0 else [],
        'mobility': fp_mz_arr[:, 1] if len(fp_mz_arr) > 0 else [],
        'snr': sn[fp_mask],
        'match_ratio': pulse_ratio[fp_mask],
        'rsd_subs': rsd_subs[fp_mask],
        'cell_occur_ratio': cell_occur_ratio[fp_mask],
        'score': score[fp_mask],
    }).to_excel(os.path.join(path_p, 'cell_specific_peak.xlsx'), index=False)

    pd.DataFrame({
        'mz': trash_mz_arr[:, 0] if len(trash_mz_arr) > 0 else [],
        'mobility': trash_mz_arr[:, 1] if len(trash_mz_arr) > 0 else [],
        'snr': sn[~fp_mask],
        'match_ratio': pulse_ratio[~fp_mask],
        'rsd_subs': rsd_subs[~fp_mask],
        'cell_occur_ratio': cell_occur_ratio[~fp_mask],
        'score': score[~fp_mask],
    }).to_excel(os.path.join(path_p, 'cell_nonspecific_peak.xlsx'), index=False)

    all_mz_arr = np.array([list(x) for x in mz_tar])
    pd.DataFrame({
        'mz': all_mz_arr[:, 0],
        'mobility': all_mz_arr[:, 1],
        'snr': sn,
        'match_ratio': pulse_ratio,
        'rsd_subs': rsd_subs,
        'cell_occur_ratio': cell_occur_ratio,
        'score': score,
    }).to_excel(os.path.join(path_p, 'all_candidate_peak.xlsx'), index=False)

    # 绘图（[FIX] 传入 ppm_frac 而非 ppm）
    plot_eics_tims(path_p, np.array(RT), eic_matrix_p, np.array(mz_tar, dtype=object),
                   (marker_mz, marker_ccs), ppm_frac, ccs_range,
                   fp_mz, fp_sn, trash_mz, trash_score,
                   plot_option, plot_rt_min, plot_rt_max)
    print("[12/14] Finish", flush=True)
    print(f"MPCell_* tables  : {output_dir}", flush=True)
    # --------------------------------------------------------
    # [13/14] Cell-specific peak 定量匹配（原 tims_quantification.ipynb）
    #   输入：cell_specific_peak（内存 fp_mz）
    #         01_feature_table.csv / 02_isotope_annotation_table.csv / 03_metabolite_annotation_table.csv
    #   输出：MPCell_01_feature_table.csv / MPCell_02_isotope_annotation_table.csv / MPCell_03_metabolite_annotation_table.csv
    # --------------------------------------------------------
    print("[13/14] Match cell-specific peaks to annotation tables", flush=True)

    # 优先使用内存里 [12/14] 计算出的 fp_mz
    cell_specific_peak = list(fp_mz)
    if len(cell_specific_peak) == 0:
        csp_path = os.path.join(output_dir, 'cell_specific_peak.npy')
        if os.path.exists(csp_path):
            cell_specific_peak = list(np.load(csp_path, allow_pickle=True))
    print(f"     cell_specific_peak count = {len(cell_specific_peak)}", flush=True)

    def _load_table(fname):
        """依次尝试 input_dir / output_dir 找表"""
        for p in (os.path.join(input_dir, fname),
                  os.path.join(output_dir, fname)):
            if os.path.exists(p):
                try:
                    return pd.read_csv(p), p
                except Exception as e:
                    print(f"     Warning: failed to read {p}: {e}", flush=True)
        return None, None

    def _match_table(fname, tag):
        df_in, src = _load_table(fname)
        if df_in is None:
            print(f"     Warning: {fname} not found, skip", flush=True)
            return
        if 'mz' not in df_in.columns or 'mobility' not in df_in.columns:
            print(f"     Warning: {fname} missing 'mz' / 'mobility' column, skip", flush=True)
            return
        mz_vals  = df_in['mz'].values
        mob_vals = df_in['mobility'].values
        idx_list = []
        for i in range(len(df_in)):
            peak_mz  = float(mz_vals[i])
            peak_mob = float(mob_vals[i])
            for peak0 in cell_specific_peak:
                mz0, ccs0 = float(peak0[0]), float(peak0[1])
                if abs(peak_mz - mz0) / mz0 < ppm_frac and abs(peak_mob - ccs0) < ccs_range:
                    idx_list.append(i)
                    break
        df_out = df_in.iloc[idx_list, :].copy()
        out_path = os.path.join(output_dir, f'MPCell_{tag}')
        df_out.to_csv(out_path, index=False)
        print(f"     [{tag}] matched {len(idx_list)} / {len(df_in)}  (src: {src})", flush=True)
        print(f"             -> {out_path}", flush=True)

    # 注：原 notebook 里 02/03 分支误用了 df1，这里分别使用正确的表
    _match_table('01_feature_table.csv',            '01_feature_table.csv')
    _match_table('02_isotope_annotation_table.csv', '02_isotope_annotation_table.csv')
    _match_table('03_metabolite_annotation_table.csv', '03_metabolite_annotation_table.csv')

    print("[13/14] Finish", flush=True)
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