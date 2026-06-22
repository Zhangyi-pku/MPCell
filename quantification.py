#!/usr/bin/env python
# coding: utf-8

# In[3]:


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
    # 归一化（以单同位素为 1）
    return np.array(intensities)

# ------------------------------------------------------------
# 主函数：同位素指认（兼容 tims / no tims 模式）
# ------------------------------------------------------------
def annotate_isotope(df,
                     mz_tol_ppm=20e-6,
                     ccs_tol=0.01,
                     isotope_delta=1.003355,
                     isotope_max_num=4,
                     isotope_int_ratio_check=True,
                     isotope_int_ratio_cutoff=500.0,
                     res_define_at=200,
                     use_ccs=None):
    df = df.copy()
    has_ccs = (use_ccs if use_ccs is not None else ('ccs' in df.columns))
    df['base_peak'] = pd.NA
    df['isotope'] = pd.NA
    df = df.sort_values('mean_intensity', ascending=False).reset_index(drop=True)
    df.insert(0, 'name', df['mz'].tolist())
    
    unmarked = df['base_peak'].isna()
    total = unmarked.sum()
    processed = 0
    print(f"Starting isotope annotation, {total} candidate peak in total", flush=True)
    
    for i, row in df[unmarked].iterrows():
        # 进度打印
        processed += 1
        if processed % max(1, total // 20) == 0 or processed == total:
            percent = processed / total * 100
            print(f"Schedule of isotope annotation: {processed}/{total} ({percent:.1f}%)", flush=True)
        
        if not pd.isna(row['base_peak']):
            continue
        base_name = row['name']
        base_mz = row['mz']
        base_intensity = row['mean_intensity']
        base_ccs = row.get('ccs') if has_ccs else None
        
        mz_isotope_list = [base_mz + i * isotope_delta for i in range(isotope_max_num)]
        mz_windows = []
        for iso_mz in mz_isotope_list:
            tol = ppm2dalton(iso_mz, mz_tol_ppm, res_define_at)
            mz_windows.append((iso_mz - tol, iso_mz + tol))
        
        candidates = []
        for j, cand in df[unmarked].iterrows():
            if cand['name'] == base_name:
                continue
            for k, (low, high) in enumerate(mz_windows):
                if low <= cand['mz'] <= high:
                    if has_ccs:
                        if 'ccs' not in cand:
                            continue
                        if abs(cand['ccs'] - base_ccs) > ccs_tol:
                            continue
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
    if 'ppm' in params:
        ppm = params.get('ppm', 10e-6)
        print(f"(Read from paras_file) ppm = {params.get('ppm')}", flush=True)
    else:
        ppm = params.get('mzppm', 20e-6)
        print(f"(Read from paras_file) ppm = {params.get('mzppm')}", flush=True)
    ccsppm = params.get('ccsppm', 0.05)               # tims 模式下的 CCS 容差
    data_type = params.get('data_type', 'auto')       # 'auto', 'no_tims', 'tims'

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
    fp_excel = os.path.join(path_p, 'c_peak_fp.xlsx')   # 仅用于参照列，但不再使用

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
    # 1. 读取用户提供的参照峰列表（Excel 文件）
    # ------------------------------------------------------------
    try:
        df_user = pd.read_excel(mz_list_path)
        # 确定模式：优先使用传入的 data_type，否则自动检测列
        if data_type == 'tims':
            # tims 模式：需要 'ccs' 和 'm/z' 两列
            if 'ccs' not in df_user.columns or ('m/z' not in df_user.columns and 'mz' not in df_user.columns):
                print("Error! Excel must include 'ccs' and 'm/z' columns under IM-based mode", file=sys.stderr, flush=True)
                sys.exit(1)
            mz_col = 'm/z' if 'm/z' in df_user.columns else 'mz'
            mz_list = df_user[mz_col].tolist()
            ccs_list = df_user['ccs'].tolist()
            ref_peaks = list(zip(ccs_list, mz_list))
            has_ccs = True
            print(f"IM-based mode：read {len(ref_peaks)} peaks (CCS, m/z)", flush=True)
        elif data_type == 'no_tims':
            # no tims 模式：只需要 m/z 列
            if 'm/z' not in df_user.columns and 'mz' not in df_user.columns:
                print("Error! Excel must include 'm/z' or 'mz' columns under no_IM mode", file=sys.stderr, flush=True)
                sys.exit(1)
            mz_col = 'm/z' if 'm/z' in df_user.columns else 'mz'
            mz_list = df_user[mz_col].tolist()
            ref_peaks = [(None, mz) for mz in mz_list]
            has_ccs = False
            print(f"No_IM mode：read {len(ref_peaks)} peaks (m/z)", flush=True)
        else:  # auto 模式：自动检测
            has_ccs = 'ccs' in df_user.columns
            mz_col = 'm/z' if 'm/z' in df_user.columns else 'mz'
            if mz_col not in df_user.columns:
                print("Error! 'm/z' column missing in excel", file=sys.stderr, flush=True)
                sys.exit(1)
            mz_list = df_user[mz_col].tolist()
            if has_ccs:
                ccs_list = df_user['ccs'].tolist()
                ref_peaks = list(zip(ccs_list, mz_list))
                print(f"Automatically detection, CCS column detected, read {len(ref_peaks)} peaks (CCS, m/z)", flush=True)
            else:
                ref_peaks = [(None, mz) for mz in mz_list]
                print(f"Automatically detection, no CCS column detected，read {len(ref_peaks)} peaks (m/z)", flush=True)
    except Exception as e:
        print(f"Failed to read provided file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # ------------------------------------------------------------
    # 2. 读取 c_peak.npy（已处理好的峰列表）
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

    for idx_ref, (ccs_ref, mz_ref) in enumerate(ref_peaks):
        mz_tol = mz_ref * mzppm
        found = False
        for idx_peak, peak_info in enumerate(c_peak):
            # 解析 c_peak 结构
            # no tims: [mz, rt_list, int_list]
            # tims:    [[ccs, mz], rt_list, int_list]
            if isinstance(peak_info[0], (list, np.ndarray)):
                # tims 模式
                ccs_peak, mz_peak = peak_info[0]
            else:
                # no tims 模式
                mz_peak = peak_info[0]
                ccs_peak = None

            # m/z 匹配
            if abs(mz_ref - mz_peak) <= mz_tol:
                # 如果参照峰有 CCS 且当前峰也有 CCS，则检查 CCS 匹配
                if has_ccs and ccs_ref is not None and ccs_peak is not None:
                    if abs(ccs_ref - ccs_peak) > ccsppm:
                        continue   # CCS 不匹配，跳过
                # 匹配成功
                intensity_list = peak_info[2]   # 原始强度列表
                blank_val = blank_aver_intensity[idx_peak]
                peak_intensity = np.array(intensity_list) - blank_val
                peak_intensity[peak_intensity < 0] = 0
                intensity_matrix.append(peak_intensity.tolist())
                found = True
                matched_count += 1
                break
        if not found:
            print(f"Error! No matched peak m/z = {mz_ref}" + (f", ccs = {ccs_ref}" if has_ccs and ccs_ref else ""), flush=True)
            unmatched.append(mz_ref)
            intensity_matrix.append([])   # 空列表占位

    print(f" {matched_count} / {len(ref_peaks)} peak matched successfully", flush=True)
   
    # ------------------------------------------------------------
    # 5. 保存强度矩阵 Excel
    # ------------------------------------------------------------
    # 行名使用 m/z（若有 CCS 可显示组合，这里保持仅 m/z）
    df_intensity = pd.DataFrame(intensity_matrix)
    mz_values = [f"{ref[1]}" for ref in ref_peaks]
    mz_values=[float(x) for x in mz_values]
    df_intensity.insert(0, 'mz', mz_values)
    ccs_values = [ref[0] for ref in ref_peaks]
    if any(x is not None for x in ccs_values):
        df_intensity.insert(0, 'ccs', ccs_values)   
    output_excel = os.path.join(output_dir, 'intensity_matrix.xlsx')
    try:
        df_intensity.to_excel(output_excel)
        print(f"Intenisity matrix saved in: {output_excel}", flush=True)
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
    ref_mean_intensity=[]
    intensity_matrix=np.array(intensity_matrix)
    for idx, intensity_list in enumerate(intensity_matrix):
        intensity_nonzero=intensity_list[intensity_list != 0]
        if len(intensity_nonzero) != 0:
            ref_mean_intensity.append(np.mean(intensity_nonzero))
        else:
            ref_mean_intensity.append(0)
    df_intensity['mean_intensity']=ref_mean_intensity
    
    if data_type=='tims':
        iso_result = annotate_isotope(df_intensity,mz_tol_ppm=mzppm,ccs_tol=0.01,use_ccs=True)
    elif data_type=='no tims':
        iso_result = annotate_isotope(df_intensity,mz_tol_ppm=mzppm,ccs_tol=0.01,use_ccs=False)
    else:
        iso_result = annotate_isotope(df_intensity,mz_tol_ppm=mzppm,ccs_tol=0.01)
    iso_output_excel = os.path.join(output_dir, 'intensity_matrix_iso_anno.xlsx')
    try:
        iso_result.to_excel(iso_output_excel)
        print(f"Intensit maxtrix with isotope annotation saved in: {iso_output_excel}", flush=True)
    except Exception as e:
        print(f"Failed to save excel file: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    print("Quantification finished", flush=True)

if __name__ == "__main__":
    main()


# In[ ]:




