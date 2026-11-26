#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8
import sys
import os
import json
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')


def main():
    if len(sys.argv) < 2:
        print("Error: no params file path given", file=sys.stderr, flush=True)
        sys.exit(1)

    with open(sys.argv[1], 'r', encoding='utf-8') as f:
        params = json.load(f)

    output_dir = params.get('output_dir', '')
    if not output_dir or not os.path.isdir(output_dir):
        print(f"Error: invalid output_dir: {output_dir}", file=sys.stderr, flush=True)
        sys.exit(1)

    print("\n========== Starting tims quantification ==========", flush=True)
    print(f"output_dir = {output_dir}", flush=True)

    cell_specific_peak_path = os.path.join(output_dir, 'cell_specific_peak.npy')
    df1_path = os.path.join(output_dir, '01_feature_table.csv')
    df2_path = os.path.join(output_dir, '02_isotope_annotation_table.csv')
    df3_path = os.path.join(output_dir, '03_metabolite_annotation_table.csv')

    for p in (cell_specific_peak_path, df1_path, df2_path, df3_path):
        if not os.path.exists(p):
            print(f"Error: required file not found: {p}", file=sys.stderr, flush=True)
            sys.exit(1)

    cell_specific_peak = np.load(cell_specific_peak_path, allow_pickle=True)
    print(f"  cell_specific_peak count: {len(cell_specific_peak)}", flush=True)

    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)
    df3 = pd.read_csv(df3_path)

    ppm_frac = 20e-6
    mob_tol = 0.05

    def _match(df, tag):
        if 'mz' not in df.columns or 'mobility' not in df.columns:
            print(f"Warning: {tag} missing 'mz' or 'mobility' column, skip", flush=True)
            return
        mz = df['mz'].values
        mob = df['mobility'].values
        idx_list = []
        for idx in range(len(df)):
            peak_mz = float(mz[idx])
            peak_mob = float(mob[idx])
            for peak0 in cell_specific_peak:
                if abs(peak_mz - peak0[0]) / peak0[0] < ppm_frac and abs(peak_mob - peak0[1]) < mob_tol:
                    idx_list.append(idx)
                    break
        df_new = df.iloc[idx_list, :].copy()
        out_path = os.path.join(output_dir, f'MPCell_{tag}')
        df_new.to_csv(out_path, index=False)
        print(f"[{tag}] matched {len(idx_list)} / {len(df)} -> {out_path}", flush=True)

    # 注意：notebook 中 02/03 误用了 df1，这里修正为 df2/df3
    _match(df1, '01_feature_table.csv')
    _match(df2, '02_isotope_annotation_table.csv')
    _match(df3, '03_metabolite_annotation_table.csv')

    print("\n*** tims quantification completed successfully ***\n", flush=True)


if __name__ == '__main__':
    main()

