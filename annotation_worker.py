#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import pandas as pd
import zlib

# 强制实时输出
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
else:
    import functools
    print = functools.partial(print, flush=True)

def load_embedded_database(name, ion_mode):
    """从内嵌的二进制文件加载数据库，返回 DataFrame，并将加合物列重命名为 'mz_加合物' 格式"""
    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS
    else:
        base_path = os.path.dirname(__file__)
    
    sheet_key = 'Positive' if ion_mode == 'Positive' else 'Negative'
    bin_filename = f"db_{name}_{sheet_key}.bin"
    bin_path = os.path.join(base_path, 'db_data_bin', bin_filename)
    
    if not os.path.exists(bin_path):
        return None
    try:
        with open(bin_path, 'rb') as f:
            compressed = f.read()
        json_str = zlib.decompress(compressed).decode('utf-8')
        records = json.loads(json_str)
        df = pd.DataFrame(records)
        # 将加合物列（以 '[' 开头）重命名为 'mz_原列名'
        rename_dict = {col: f'mz_{col}' for col in df.columns if col.startswith('[')}
        df.rename(columns=rename_dict, inplace=True)
        return df
    except Exception as e:
        print(f"Failed to load embedded database {name}: {e}", file=sys.stderr, flush=True)
        return None

def ppm_to_abs(mz, ppm):
    return mz * ppm * 1e-6

def match_single_mz_all(mz, db_df, adduct_priority, ppm_tol):
    matches = []
    for adduct in adduct_priority:
        col_name = f'mz_{adduct}'
        if col_name not in db_df.columns:
            continue
        tol = ppm_to_abs(mz, ppm_tol)
        low, high = mz - tol, mz + tol
        mask = (db_df[col_name] >= low) & (db_df[col_name] <= high)
        for _, row in db_df[mask].iterrows():
            theo_mz = row[col_name]
            ppm_error = abs((theo_mz - mz) / mz) * 1e6
            matches.append({
                'adduct': adduct,
                'theoretical_mz': theo_mz,
                'name': row.get('Name', None),
                'formula': row.get('Formula', None),
                'hmdb_id': row.get('HMDB', None),
                'ppm_error': ppm_error
            })
    priority_map = {ad: i for i, ad in enumerate(adduct_priority)}
    matches.sort(key=lambda x: (priority_map.get(x['adduct'], len(adduct_priority)), x['ppm_error']))
    return matches

def main():
    if len(sys.argv) < 2:
        print("Error! No paras_file path specified", file=sys.stderr, flush=True)
        sys.exit(1)

    param_file = sys.argv[1]
    with open(param_file, 'r', encoding='utf-8') as f:
        params = json.load(f)

    input_mz_file = params.get('input_mz_file')
    output_dir = params.get('output_dir')
    db_files = params.get('db_files', [])          # 列表，每个元素为 {'name': ..., 'path': ...}
    adduct_priority = params.get('adduct_priority', [])
    ppm_tol = params.get('ppm_tol', 10)
    max_candidates = params.get('max_candidates', 5)
    ion_mode = params.get('ion_mode', 'Positive')
    db_priority = params.get('db_priority', [])   # 数据库名称优先级列表

    # 根据离子模式过滤加合物
    allowed = {
        '[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+','[M+H-H2O]+','[M+H-2H2O]+','[M+2Na-H]+','[M+2K-H]+','[M+2Na]2+',
        '[M+2K]2+','[M+2H]2+','[M+3H]3+','[M+4H]4+','[M.]+'
    } if ion_mode == 'Positive' else {
        '[M-H]-','[M+Cl]-','[M+Na-2H]-','[M+OAc]-','[M+HCOO]-','[M+HCO3]-','[M-2H]2-',
        '[M-3H]3-','[M-4H]4-','[M-CH3]-','[M.]-'
    }
    filtered_priority = [ad for ad in adduct_priority if ad in allowed]
    if not filtered_priority:
        print(f"Warning! No ion adducts selected in the {ion_mode} ion mode, default list used", flush=True)
        filtered_priority = list(allowed)
    print(f"Ion poalrity mode: {ion_mode}, ion adduct priority: {filtered_priority}", flush=True)

    # 读取检测 m/z 文件
    try:
        df_mz = pd.read_excel(input_mz_file)
        if 'isotope' in df_mz.columns:
            df_mz = df_mz[df_mz['isotope'] == '[M]'].copy()
            print(f"{len(df_mz)} M0 peaks reserved after isotope filtration", flush=True)
        if 'mz' in df_mz.columns:
            mz_list = df_mz['mz'].tolist()
        elif 'm/z' in df_mz.columns:
            mz_list = df_mz['m/z'].tolist()
        else:
            mz_list = df_mz.iloc[:, 0].tolist()
        print(f"Read {len(mz_list)} m/z values", flush=True)
    except Exception as e:
        print(f"Failed to read input files: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

    # 加载数据库
    embedded_names = ["lcms_mda", "lcms_hela", "lcms_ninecells", "sc_paper", "hmdb", "lipidmaps"]
    db_dict = {}

    # 先加载内嵌数据库
    for name in embedded_names:
        df = load_embedded_database(name, ion_mode)
        if df is not None:
            db_dict[name] = df
            print(f"Emebedded database {name} ({ion_mode} ion mode) successfully loaded, {len(df)} record", flush=True)

    # 再加载自定义数据库（如果提供）
    for db_cfg in db_files:
        db_name = db_cfg.get('name')
        db_path = db_cfg.get('path')
        # 跳过已加载的内嵌名称（避免重复）
        if db_name in db_dict:
            continue
        if not db_path or not os.path.exists(db_path):
            print(f"Warning! Self-defined database does not existed in {db_path}", flush=True)
            continue
        try:
            # 尝试按离子模式读取工作表
            df_custom = pd.read_excel(db_path, sheet_name=ion_mode)
            # 重命名加合物列
            rename_dict = {col: f'mz_{col}' for col in df_custom.columns if col.startswith('[')}
            df_custom.rename(columns=rename_dict, inplace=True)
            db_dict[db_name] = df_custom
            print(f"Self-defined database {db_name} ({ion_mode} ion mode) successfully loaded, {len(df_custom)} record", flush=True)
        except Exception:
            # 若没有对应工作表，则读取整个文件
            try:
                df_custom = pd.read_excel(db_path)
                rename_dict = {col: f'mz_{col}' for col in df_custom.columns if col.startswith('[')}
                df_custom.rename(columns=rename_dict, inplace=True)
                db_dict[db_name] = df_custom
                print(f"Self-defined database {db_name} successfully loaded (defual sheet), {len(df_custom)} record", flush=True)
            except Exception as e2:
                print(f"Failed to read self-defined database{db_name}: {e2}", flush=True)

    if not db_dict:
        print("Error! No database successfully loaded", file=sys.stderr, flush=True)
        sys.exit(1)

    # 确定数据库优先级顺序
    if db_priority:
        final_priority = [db for db in db_priority if db in db_dict]
        for db in db_dict:
            if db not in final_priority:
                final_priority.append(db)
    else:
        final_priority = list(db_dict.keys())
    print(f"Database priority: {final_priority}", flush=True)

    # 匹配
    all_results = []
    total = len(mz_list)
    for idx, mz in enumerate(mz_list):
        if idx % 100 == 0:
            print(f"Schedule: {idx}/{total}", flush=True)
        found = False
        best_candidates = []
        for db_name in final_priority:
            db_df = db_dict[db_name]
            matches = match_single_mz_all(mz, db_df, filtered_priority, ppm_tol)
            if matches:
                for m in matches:
                    m['database'] = db_name
                best_candidates = matches[:max_candidates]
                found = True
                break
        if not found:
            all_results.append({
                'Detected_mz': mz,
                'Database': None,
                'Adduct': None,
                'Theoretical_mz': None,
                'Name': None,
                'Formula': None,
                'HMDB_ID': None,
                'ppm_error': None,
                'Rank': 1
            })
        else:
            for rank, cand in enumerate(best_candidates, start=1):
                all_results.append({
                    'Detected_mz': mz,
                    'Database': cand['database'],
                    'Adduct': cand['adduct'],
                    'Theoretical_mz': cand['theoretical_mz'],
                    'Name': cand['name'],
                    'Formula': cand['formula'],
                    'HMDB_ID': cand['hmdb_id'],
                    'ppm_error': cand['ppm_error'],
                    'Rank': rank
                })

    # 保存结果
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'annotation_results.xlsx')
    df_result = pd.DataFrame(all_results)
    df_result.to_excel(output_file, index=False)
    print(f"Annotation results saved in: {output_file}", flush=True)
    print("Finish", flush=True)

if __name__ == "__main__":
    main()

