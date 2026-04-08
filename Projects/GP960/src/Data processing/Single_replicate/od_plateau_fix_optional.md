```python
import sys
import os
import pandas as pd
import numpy as np

def detect_first_peak(y, drop_len=3):
    y = np.asarray(y, dtype=float)
    n = len(y)
    run = 0
    for i in range(1, n):
        a = y[i-1]
        b = y[i]
        if np.isnan(a) or np.isnan(b):
            run = 0
            continue
        if b < a:
            run += 1
        else:
            run = 0
        if run >= drop_len:
            idx = i - drop_len
            return idx, float(y[idx])
    return int(np.nanargmax(y)), float(np.nanmax(y))

def fix_series(y, peak_idx):
    y = np.asarray(y, dtype=float)
    yf = y.copy()
    if peak_idx is None:
        return yf
    val = yf[peak_idx]
    yf[peak_idx:] = val
    return yf

def process_csv(input_csv, output_dir, drop_len=3):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(input_csv, index_col=0)
    results = []
    fixed = {}
    for col in df.columns:
        y = df[col].values
        peak_idx, peak_y = detect_first_peak(y, drop_len=drop_len)
        peak_time = df.index[peak_idx]
        fixed[col] = fix_series(y, peak_idx)
        results.append({'gene': col, 'peak_time': peak_time, 'peak_y': peak_y})
    coord_df = pd.DataFrame(results)
    #cp = os.path.join(output_dir, 'od_first_peak_coords.csv')
    #coord_df.to_csv(cp, index=False)
    print("peak coords:", coord_df)
    fixed_df = pd.DataFrame(fixed, index=df.index)
    file_name = input_csv.rsplit('/', 1)[-1]
    fp = os.path.join(output_dir, file_name)
    fixed_df.to_csv(fp)
    return fp

if __name__ == '__main__':
    ### MSG
    output_dir = '/data/zuoll/1.project/02.GOE/02_result/02.MSG/01.1.ppraw_data_od_fix'
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/02.MSG/01.ppraw_data/01.Goe_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/02.MSG/01.ppraw_data/02.Goe_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/02.MSG/01.ppraw_data/03.Ctrol_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/02.MSG/01.ppraw_data/04.Ctrol_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    ### NaCL
    output_dir = '/data/zuoll/1.project/02.GOE/02_result/03.NaCl/01.1.ppraw_data_od_fix'
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/03.NaCl/01.ppraw_data/01.Goe_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/03.NaCl/01.ppraw_data/02.Goe_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/03.NaCl/01.ppraw_data/03.Ctrol_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/03.NaCl/01.ppraw_data/04.Ctrol_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    ### AmGlu 97tp
    output_dir = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.1.ppraw_data_od_fix/97tp'
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/97tp/01.Goe_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/97tp/02.Goe_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/97tp/03.Ctrol_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/97tp/04.Ctrol_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    ### AmGlu 145tp
    output_dir = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.1.ppraw_data_od_fix/145tp'
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/145tp/01.Goe_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/145tp/02.Goe_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/145tp/03.Ctrol_stress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)
    input_csv = '/data/zuoll/1.project/02.GOE/02_result/01.AmGlu/01.ppraw_data/145tp/04.Ctrol_Nonstress_OD.csv'
    fp = process_csv(input_csv, output_dir, drop_len=5)

```