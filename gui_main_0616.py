#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import subprocess
import sys
import json
import threading
import os
import queue

# ------------------------------------------------------------
# 全局变量
# ------------------------------------------------------------
current_proc = None
plot_var = None
filter_mode_var = None
entry_threshold = None
threshold_hint = None
threshold_is_default = True
mode_var = None
ion_mode_var = None
orig_rt_unit_var = None
quant_proc = None  
anno_proc = None   
embedded_db_vars = {}   # 存储内嵌数据库勾选状态
db_priority_list = []   # 存储已启用的数据库名称列表（按优先级顺序）
custom_db_path = None   # 存储自定义数据库文件路径
adduct_priority_list = []   # 存储已启用的加合物列表（按优先级顺序）

# 背景标记动态行列表
bg_rows = []

# ------------------------------------------------------------
# 文件/文件夹浏览回调函数
# ------------------------------------------------------------
def browse_input_file(entry_widget):
    path = filedialog.askopenfilename(
        title="Select Input File",
        filetypes=[("mzML files", "*.mzml"), ("All files", "*.*")]
    )
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)
        if not entry_output_dir.get().strip():
            default_output_dir = get_default_output_dir(path)
            entry_output_dir.delete(0, tk.END)
            entry_output_dir.insert(0, default_output_dir)

def browse_output_dir(entry_widget):
    path = filedialog.askdirectory(title="Select Output Folder")
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)

def browse_mz_list(entry_widget):
    path = filedialog.askopenfilename(
        title="Select peak_list Excel File",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
    )
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)

def browse_mz_file(entry_widget):
    path = filedialog.askopenfilename(title="Select m/z Input File (Excel)", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)

def browse_db_file(entry_widget):
    path = filedialog.askopenfilename(title="Select Database File", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)

def update_db_priority_list():
    """根据内嵌数据库勾选状态和自定义数据库，重建优先级列表（保留原顺序）"""
    global db_priority_list, custom_db_path
    # 收集当前启用的内嵌数据库名称
    new_embedded = [name for name, var in embedded_db_vars.items() if var.get()]
    # 自定义数据库（仅当路径存在时）
    new_custom = ["Custom"] if custom_db_path else []
    # 保留原顺序中仍然存在的项
    old_order = db_priority_list.copy()
    result = []
    for db in old_order:
        if db == "Custom":
            # 仅当 custom_db_path 非空时才保留 Custom
            if custom_db_path:
                result.append(db)
        elif db in new_embedded:
            result.append(db)
    # 添加新增的内嵌数据库（之前未在结果中的）
    for db in new_embedded:
        if db not in result:
            result.append(db)
    # 添加自定义数据库（如果尚未添加）
    if new_custom and "Custom" not in result:
        result.append("Custom")
    db_priority_list = result
    refresh_db_listbox()

def refresh_db_listbox():
    """刷新优先级列表框显示"""
    listbox_db_priority.delete(0, tk.END)
    for db in db_priority_list:
        listbox_db_priority.insert(tk.END, db)

def add_custom_db():
    """添加自定义数据库文件"""
    global custom_db_path
    path = filedialog.askopenfilename(title="Select Custom Database", filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
    if path:
        custom_db_path = path
        entry_db_file.delete(0, tk.END)
        entry_db_file.insert(0, path)
        update_db_priority_list()

def remove_custom_db():
    """移除自定义数据库"""
    global custom_db_path
    custom_db_path = None
    entry_db_file.delete(0, tk.END)
    update_db_priority_list()
    listbox_db_priority.selection_clear(0, tk.END)

def add_db():
    """从内嵌数据库中选择未添加的数据库添加到优先级列表"""
    current_embedded = [name for name in embedded_names if name in db_priority_list]
    available = [name for name in embedded_names if name not in current_embedded]
    if not available:
        messagebox.showinfo("Info", "All embedded databases have been added")
        return
    win = tk.Toplevel(root)
    win.title("Select Database")
    win.geometry("300x200")
    win.transient(root)
    win.grab_set()
    lb = tk.Listbox(win, selectmode=tk.SINGLE, height=8)
    lb.pack(fill='both', expand=True, padx=10, pady=10)
    for db in available:
        lb.insert(tk.END, db)
    def confirm():
        sel = lb.curselection()
        if sel:
            selected = lb.get(sel[0])
            db_priority_list.append(selected)
            if selected in embedded_db_vars:
                embedded_db_vars[selected].set(True)
            refresh_db_listbox()
        win.destroy()
    ttk.Button(win, text="OK", command=confirm).pack(pady=5)

def remove_db():
    """从优先级列表中移除选中数据库，并更新对应复选框及自定义数据库状态"""
    sel = listbox_db_priority.curselection()
    if not sel:
        messagebox.showinfo("Info", "Please select a database to remove")
        return
    idx = sel[0]
    removed = db_priority_list[idx]
    db_priority_list.pop(idx)
    if removed in embedded_db_vars:
        embedded_db_vars[removed].set(False)
    if removed == "Custom":
        global custom_db_path
        custom_db_path = None
        entry_db_file.delete(0, tk.END)
    refresh_db_listbox()
    listbox_db_priority.selection_clear(0, tk.END)
    
def move_db_up():
    """将优先级列表中的选中项上移"""
    selected = listbox_db_priority.curselection()
    if not selected:
        return
    idx = selected[0]
    if idx > 0:
        db_priority_list[idx], db_priority_list[idx-1] = db_priority_list[idx-1], db_priority_list[idx]
        refresh_db_listbox()
        listbox_db_priority.selection_set(idx-1)

def move_db_down():
    """将优先级列表中的选中项下移"""
    selected = listbox_db_priority.curselection()
    if not selected:
        return
    idx = selected[0]
    if idx < len(db_priority_list)-1:
        db_priority_list[idx], db_priority_list[idx+1] = db_priority_list[idx+1], db_priority_list[idx]
        refresh_db_listbox()
        listbox_db_priority.selection_set(idx+1)

def run_annotation_subprocess():
    global anno_proc, db_priority_list, custom_db_path
    # 获取输入
    mz_file = entry_mz_file.get().strip()
    if not mz_file:
        messagebox.showwarning("Warning", "Please select an input m/z file")
        return
    output_dir = entry_output_dir.get().strip()
    if not output_dir:
        messagebox.showwarning("Warning", "Please select an output folder first")
        return
    try:
        ppm_tol = float(entry_ppm_anno.get())
    except:
        ppm_tol = 10
    try:
        max_candidates = int(entry_max_candidates.get())
    except:
        max_candidates = 5
    adduct_priority = get_adduct_priority()
    if not adduct_priority:
        messagebox.showwarning("Warning", "Please select at least one adduct type")
        return
    ion_mode = ion_mode_var_anno.get()

    # 数据库优先级列表（已包含用户选中的内嵌数据库名称，按优先级顺序）
    db_priority = db_priority_list
    if not db_priority:
        messagebox.showwarning("Warning", "Please select at least one embedded database")
        return

    # 可选的自定义数据库文件路径
    custom_db = entry_db_file.get().strip()
    if custom_db and not os.path.exists(custom_db):
        messagebox.showwarning("Warning", f"Custom database file does not exist: {custom_db}")
        return

    # 构建参数字典
    params = {
        'input_mz_file': mz_file,
        'output_dir': output_dir,
        'ppm_tol': ppm_tol,
        'max_candidates': max_candidates,
        'adduct_priority': adduct_priority,
        'ion_mode': ion_mode,
        'db_priority': db_priority,      # 内嵌数据库名称列表
        'custom_db': custom_db           # 自定义数据库路径（可为空）
    }
    # 写入临时参数文件
    param_file_path = os.path.join(output_dir, 'temp_annotation_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Unable to create parameter file: {e}")
        return

    # 启动子进程（以下代码保持不变）
    btn_anno_run.config(state="disabled")
    btn_anno_terminate.config(state="normal")
    anno_result_text.delete(1.0, tk.END)
    anno_result_text.insert(tk.END, "Starting annotation...\n")
    root.update()

    output_queue = queue.Queue()

    def read_stream(stream, stream_name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((stream_name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((stream_name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global anno_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            cmd = get_annotation_command()
            proc = subprocess.Popen(
                cmd + [param_file_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                env=env
            )
            anno_proc = proc

            t_stdout = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_stderr = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_stdout.daemon = True
            t_stderr.daemon = True
            t_stdout.start()
            t_stderr.start()

            return_code = proc.wait()
            t_stdout.join(timeout=1)
            t_stderr.join(timeout=1)

            output_queue.put(("DONE", return_code))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    return_code = item[1]
                    if return_code == 0:
                        anno_result_text.insert(tk.END, "\n*** Annotation completed successfully ***\n")
                    else:
                        anno_result_text.insert(tk.END, f"\n*** Annotation exited with error code {return_code} ***\n")
                    try:
                        os.remove(param_file_path)
                    except:
                        pass
                    btn_anno_run.config(state="normal")
                    btn_anno_terminate.config(state="disabled")
                    return
                elif item[0] == "ERROR":
                    anno_result_text.insert(tk.END, f"\n!!! Error: {item[1]} !!!\n")
                else:
                    prefix = ">> " if item[0] == "STDOUT" else "!! "
                    anno_result_text.insert(tk.END, f"{prefix}{item[1]}\n")
                    anno_result_text.see(tk.END)
        except queue.Empty:
            pass
        root.after(100, poll_queue)

    threading.Thread(target=target, daemon=True).start()
    root.after(100, poll_queue)

def terminate_annotation_subprocess():
    global anno_proc
    if anno_proc is not None:
        try:
            anno_proc.terminate()
            try:
                anno_proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                anno_proc.kill()
            anno_result_text.insert(tk.END, "\n⚠️ User terminated annotation manually\n")
        except Exception as e:
            messagebox.showerror("Termination Error", str(e))
        finally:
            anno_proc = None
            btn_anno_run.config(state="normal")
            btn_anno_terminate.config(state="disabled")
    else:
        messagebox.showinfo("Info", "No annotation process running")

def get_default_output_dir(input_file):
    if not input_file:
        return ""
    if mode_var.get() == "tims":
        base_name = os.path.basename(input_file.rstrip(os.sep))
        dir_name = os.path.dirname(input_file)
    else:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        dir_name = os.path.dirname(input_file)
    return os.path.join(dir_name, base_name)

def get_worker_command():
    if getattr(sys, 'frozen', False):
        base_dir = os.path.dirname(sys.executable)
        if mode_var.get() == "tims":
            worker_name = "worker_tims.exe"
        else:
            worker_name = "worker.exe"
        worker_path = os.path.join(base_dir, worker_name)
        return [worker_path]
    else:
        if mode_var.get() == "tims":
            return [sys.executable, "worker_tims.py"]
        else:
            return [sys.executable, "worker.py"]

def get_quant_command():
    if getattr(sys, 'frozen', False):
        base_dir = os.path.dirname(sys.executable)
        quant_path = os.path.join(base_dir, 'quantification.exe')
        return [quant_path]
    else:
        return [sys.executable, "quantification.py"]

def get_annotation_command():
    if getattr(sys, 'frozen', False):
        base_dir = os.path.dirname(sys.executable)
        worker_path = os.path.join(base_dir, 'annotation_worker.exe')
        return [worker_path]
    else:
        return [sys.executable, "annotation_worker.py"]

# ------------------------------------------------------------
# 背景标记行管理函数（根据模式决定是否显示 CCS 列）
# ------------------------------------------------------------
def add_bg_row(ccs_val="", mz_val="", int_val=""):
    row_frame = tk.Frame(bg_rows_container)
    row_frame.pack(fill='x', pady=2)

    entries = {}
    if mode_var.get() == "tims":
        entry_ccs = tk.Entry(row_frame, width=12)
        entry_ccs.insert(0, ccs_val)
        entry_ccs.pack(side='left', padx=2)
        entries['ccs'] = entry_ccs
    else:
        entries['ccs'] = None

    entry_mz = tk.Entry(row_frame, width=15)
    entry_mz.insert(0, mz_val)
    entry_mz.pack(side='left', padx=2)

    entry_int = tk.Entry(row_frame, width=12)
    entry_int.insert(0, int_val)
    entry_int.pack(side='left', padx=2)

    # 删除按钮明确设置宽度和文字
    del_btn = tk.Button(row_frame, text="✖", command=lambda: delete_bg_row(row_frame))
    del_btn.pack(side='left', padx=5)

    bg_rows.append({
        'frame': row_frame,
        'ccs': entries.get('ccs', None),
        'mz': entry_mz,
        'intensity': entry_int
    })

def delete_bg_row(row_frame):
    for i, row in enumerate(bg_rows):
        if row['frame'] == row_frame:
            row_frame.destroy()
            bg_rows.pop(i)
            break

def clear_bg_rows():
    for row in bg_rows:
        row['frame'].destroy()
    bg_rows.clear()

def populate_bg_rows(defaults):
    clear_bg_rows()
    for item in defaults:
        if mode_var.get() == "tims":
            ccs, mz, ints = item
            add_bg_row(ccs, mz, ints)
        else:
            # no tims 模式：defaults 中每个元素为 (mz, ints)
            mz, ints = item
            add_bg_row("", mz, ints)

# ------------------------------------------------------------
# 正负模式更新参数
# ------------------------------------------------------------
def update_mode_params():
    mode = ion_mode_var.get()
    if mode_var.get() == "tims":
        if mode == "Positive":
            marker_ccs = "1.39"
            marker_mz = "760.5851"
            intensity_thres_val = "5e4"
            bg_defaults = [
                ("1.39", "760.5851", "1e4"),
                ("1.36", "762.5280", "0"),
                ("1.39", "775.5667", "0")
            ]
        else:
            marker_ccs = "1.39"
            marker_mz = "306.0767"
            intensity_thres_val = "5e4"
            bg_defaults = [
                ("1.39", "306.0767", "2e4"),
                ("1.36", "804.5732", "0"),
                ("1.39", "742.5393", "0")
            ]
        # 更新提示标签（控件显示由 on_mode_change 控制）
        label_bg_hint.config(text="Background markers (CCS, m/z, Intensity max):")
    else:  # no_tims
        if mode == "Positive":
            marker_ccs = ""
            marker_mz = "760.5851"
            intensity_thres_val = "1e6"
            bg_defaults = [
                ("760.5851", "0"),
                ("732.5532", "0")
            ]
        else:
            marker_ccs = ""
            marker_mz = "306.0767"
            intensity_thres_val = "5e6"
            bg_defaults = [
                ("306.0767", "2e6"),
                ("804.5732", "0"),
                ("742.5393", "0")
            ]
        label_bg_hint.config(text="Background markers (m/z, Intensity max):")

    entry_marker_ccs.delete(0, tk.END)
    entry_marker_ccs.insert(0, marker_ccs)
    entry_marker_mz.delete(0, tk.END)
    entry_marker_mz.insert(0, marker_mz)
    entry_intensity_thres.delete(0, tk.END)
    entry_intensity_thres.insert(0, intensity_thres_val)
    populate_bg_rows(bg_defaults)

# ------------------------------------------------------------
# 其他通用函数
# ------------------------------------------------------------
def on_threshold_edit(event):
    global threshold_is_default
    threshold_is_default = False

def update_threshold_label():
    global threshold_is_default
    mode = filter_mode_var.get()
    if mode == "snr":
        threshold_hint.config(text="(SNR threshold, default 3)")
        if threshold_is_default:
            entry_threshold.delete(0, tk.END)
            entry_threshold.insert(0, "3")
            threshold_is_default = True
    else:
        threshold_hint.config(text="(Score threshold, default 10)")
        if threshold_is_default:
            entry_threshold.delete(0, tk.END)
            entry_threshold.insert(0, "10")
            threshold_is_default = True

# ------------------------------------------------------------
# 参数收集
# ------------------------------------------------------------
def get_parameters():
    params = {}
    params['peak_intensity_thres'] = float(entry_peak_intensity.get())
    params['ppm'] = (1e-6) * float(entry_ppm.get())
    rt_min_input = float(entry_rt_min.get())
    rt_max_input = float(entry_rt_max.get())
    orig_unit = orig_rt_unit_var.get()
    if orig_unit == "second":
        params['rt_range_min'] = rt_min_input * 60.0
        params['rt_range_max'] = rt_max_input * 60.0
    else:
        params['rt_range_min'] = rt_min_input
        params['rt_range_max'] = rt_max_input
    params['cell_spec'] = int(entry_cell_spec.get())
    params['cell_spec_lim'] = int(entry_cell_spec_lim.get())
    
    # 离子模式标记
    if mode_var.get() == "tims":
        marker_ccs = float(entry_marker_ccs.get()) if entry_marker_ccs.get().strip() else 0.0
        marker_mz = float(entry_marker_mz.get())
        params['marker'] = (marker_ccs, marker_mz)
    else:
        params['marker'] = float(entry_marker_mz.get())
    params['intensity_thres'] = float(entry_intensity_thres.get())
    
    bg_marker_list = []
    bg_intensity_list = []
    for row in bg_rows:
        mz_str = row['mz'].get().strip()
        int_str = row['intensity'].get().strip()
        if mz_str == "":
            continue
        if mode_var.get() == "tims":
            ccs_str = row['ccs'].get().strip() if row['ccs'] else ""
            if ccs_str == "":
                continue
            bg_marker_list.append((float(ccs_str), float(mz_str)))
        else:
            bg_marker_list.append(float(mz_str))
        bg_intensity_list.append(float(int_str) if int_str else 0.0)
    params['bg_marker'] = bg_marker_list
    params['bg_intensity_max'] = bg_intensity_list
    
    cell_spec = params['cell_spec']
    if cell_spec % 2 == 0:
        spec_range = list(range(int(-cell_spec/2), int(cell_spec/2)))
    else:
        spec_range = list(range(int(-cell_spec/2), int(cell_spec/2) + 1))
    params['spec_range'] = spec_range
    
    params['input_file'] = entry_input_file.get().strip()
    output_dir = entry_output_dir.get().strip()
    if not output_dir and params['input_file']:
        output_dir = get_default_output_dir(params['input_file'])
    params['output_dir'] = output_dir
    params['plot_option'] = plot_var.get()
    plot_min_input = float(entry_plot_rt_min.get())
    plot_max_input = float(entry_plot_rt_max.get())
    if orig_unit == "second":
        params['plot_rt_min'] = plot_min_input * 60.0
        params['plot_rt_max'] = plot_max_input * 60.0
    else:
        params['plot_rt_min'] = plot_min_input
        params['plot_rt_max'] = plot_max_input
    params['filter_mode'] = filter_mode_var.get()
    params['filter_threshold'] = float(entry_threshold.get())
    params['data_type'] = mode_var.get()
    return params

# ------------------------------------------------------------
# 终止子进程
# ------------------------------------------------------------
def terminate_subprocess():
    global current_proc
    if current_proc is not None:
        try:
            current_proc.terminate()
            try:
                current_proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                current_proc.kill()
            result_text.insert(tk.END, "\n⚠️ User terminated subprocess manually\n")
            status_label.config(text="Terminated")
        except Exception as e:
            messagebox.showerror("Termination Error", str(e))
        finally:
            current_proc = None
            btn_run.config(state="normal")
            btn_terminate.config(state="disabled")
    else:
        messagebox.showinfo("Info", "No subprocess running")

def terminate_quant_subprocess():
    global quant_proc
    if quant_proc is not None:
        try:
            quant_proc.terminate()
            try:
                quant_proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                quant_proc.kill()
            quant_result_text.insert(tk.END, "\n⚠️ User terminated quantification manually\n")
        except Exception as e:
            messagebox.showerror("Termination Error", str(e))
        finally:
            quant_proc = None
            btn_quant_run.config(state="normal")
            btn_quant_terminate.config(state="disabled")
    else:
        messagebox.showinfo("Info", "No quantification process running")

# ------------------------------------------------------------
# 子进程运行
# ------------------------------------------------------------
def run_subprocess_threaded():
    global current_proc
    try:
        params = get_parameters()
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Failed to parse parameters: {e}")
        return

    if not params['input_file']:
        messagebox.showwarning("Warning", "Please select an input file")
        return
    if not params['output_dir']:
        messagebox.showwarning("Warning", "Please select an output folder")
        return

    os.makedirs(params['output_dir'], exist_ok=True)

    btn_run.config(state="disabled")
    btn_terminate.config(state="normal")
    status_label.config(text="Subprocess running...")
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, "Starting subprocess...\n")
    root.update()

    param_file_path = os.path.join(params['output_dir'], 'temp_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Unable to create parameter file: {e}")
        btn_run.config(state="normal")
        btn_terminate.config(state="disabled")
        return

    output_queue = queue.Queue()

    def read_stream(stream, stream_name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((stream_name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((stream_name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global current_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            worker_cmd = get_worker_command()
            proc = subprocess.Popen(
                worker_cmd + [param_file_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                env=env
            )
            current_proc = proc

            t_stdout = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_stderr = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_stdout.daemon = True
            t_stderr.daemon = True
            t_stdout.start()
            t_stderr.start()

            return_code = proc.wait()
            t_stdout.join(timeout=1)
            t_stderr.join(timeout=1)

            output_queue.put(("DONE", return_code))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    return_code = item[1]
                    if return_code == 0:
                        result_text.insert(tk.END, "\n*** Subprocess completed successfully ***\n")
                        status_label.config(text="Completed")
                    else:
                        result_text.insert(tk.END, f"\n*** Subprocess exited with error code {return_code} ***\n")
                        status_label.config(text="Error")
                    try:
                        os.remove(param_file_path)
                    except:
                        pass
                    btn_run.config(state="normal")
                    btn_terminate.config(state="disabled")
                    status_label.config(text="Ready")
                    return
                elif item[0] == "ERROR":
                    result_text.insert(tk.END, f"\n!!! Error: {item[1]} !!!\n")
                else:
                    prefix = ">> " if item[0] == "STDOUT" else "!! "
                    result_text.insert(tk.END, f"{prefix}{item[1]}\n")
                    result_text.see(tk.END)
        except queue.Empty:
            pass
        root.after(100, poll_queue)

    threading.Thread(target=target, daemon=True).start()
    root.after(100, poll_queue)

def run_quant_subprocess():
    global quant_proc
    # 获取用户输入
    mz_list_path = entry_mz_list.get().strip()
    if not mz_list_path:
        messagebox.showwarning("Warning", "Please select an mz_list file")
        return
    input_file = entry_input_file.get().strip()
    if not input_file:
        messagebox.showwarning("Warning", "Please select an input file first")
        return
    output_dir = entry_output_dir.get().strip()
    if not output_dir:
        messagebox.showwarning("Warning", "Please select an output folder first")
        return

    # 构建参数字典（注意：quantification.py 需要接收的参数可能还包括 ppm 等，请根据您的实际情况调整）
    params = {
        'mz_list': mz_list_path,
        'input_file': input_file,
        'output_dir': output_dir
    }
    # 可选：如果需要传递 ppm 或其他参数，可以添加
    try:
        ppm_val = float(entry_ppm.get()) * 1e-6
        params['ppm'] = ppm_val
    except:
        pass

    # 写入临时参数文件
    param_file_path = os.path.join(output_dir, 'temp_quant_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Unable to create parameter file: {e}")
        return

    # 禁用运行按钮，启用终止按钮
    btn_quant_run.config(state="disabled")
    btn_quant_terminate.config(state="normal")
    quant_result_text.delete(1.0, tk.END)
    quant_result_text.insert(tk.END, "Starting quantification...\n")
    root.update()

    output_queue = queue.Queue()

    def read_stream(stream, stream_name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((stream_name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((stream_name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global quant_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            cmd = get_quant_command()
            proc = subprocess.Popen(
                cmd + [param_file_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                env=env
            )
            quant_proc = proc

            t_stdout = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_stderr = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_stdout.daemon = True
            t_stderr.daemon = True
            t_stdout.start()
            t_stderr.start()

            return_code = proc.wait()
            t_stdout.join(timeout=1)
            t_stderr.join(timeout=1)

            output_queue.put(("DONE", return_code))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    return_code = item[1]
                    if return_code == 0:
                        quant_result_text.insert(tk.END, "\n*** Quantification completed successfully ***\n")
                    else:
                        quant_result_text.insert(tk.END, f"\n*** Quantification exited with error code {return_code} ***\n")
                    try:
                        os.remove(param_file_path)
                    except:
                        pass
                    btn_quant_run.config(state="normal")
                    btn_quant_terminate.config(state="disabled")
                    return
                elif item[0] == "ERROR":
                    quant_result_text.insert(tk.END, f"\n!!! Error: {item[1]} !!!\n")
                else:
                    prefix = ">> " if item[0] == "STDOUT" else "!! "
                    quant_result_text.insert(tk.END, f"{prefix}{item[1]}\n")
                    quant_result_text.see(tk.END)
        except queue.Empty:
            pass
        root.after(100, poll_queue)

    threading.Thread(target=target, daemon=True).start()
    root.after(100, poll_queue)

# ------------------------------------------------------------
# 模式切换回调
# ------------------------------------------------------------
def on_mode_change():
    mode = mode_var.get()
    if mode == "tims":
        # 显示 tims 控件
        label_ccsppm.grid()
        entry_ccsppm.grid()
        label_eim_span.grid()
        entry_eim_span.grid()
        label_ppm.config(text="mz tolerance [ppm]:")
        # 重置默认值
        entry_peak_intensity.delete(0, tk.END)
        entry_peak_intensity.insert(0, "0")
        entry_ppm.delete(0, tk.END)
        entry_ppm.insert(0, "20")
        entry_ccsppm.delete(0, tk.END)
        entry_ccsppm.insert(0, "0.05")
        entry_rt_min.delete(0, tk.END)
        entry_rt_min.insert(0, "0")
        entry_rt_max.delete(0, tk.END)
        entry_rt_max.insert(0, "15")
        entry_cell_spec.delete(0, tk.END)
        entry_cell_spec.insert(0, "1")
        entry_cell_spec_lim.delete(0, tk.END)
        entry_cell_spec_lim.insert(0, "1")
        entry_eim_span.delete(0, tk.END)
        entry_eim_span.insert(0, "27")
        # 显示 marker CCS 输入框
        label_marker_ccs.grid()
        entry_marker_ccs.grid()
    else:
        # 隐藏 tims 控件
        label_ccsppm.grid_remove()
        entry_ccsppm.grid_remove()
        label_eim_span.grid_remove()
        entry_eim_span.grid_remove()
        label_ppm.config(text="mz tolerance [ppm]:")
        # 重置默认值
        entry_peak_intensity.delete(0, tk.END)
        entry_peak_intensity.insert(0, "1e4")
        entry_ppm.delete(0, tk.END)
        entry_ppm.insert(0, "10")
        entry_rt_min.delete(0, tk.END)
        entry_rt_min.insert(0, "0")
        entry_rt_max.delete(0, tk.END)
        entry_rt_max.insert(0, "100")
        entry_cell_spec.delete(0, tk.END)
        entry_cell_spec.insert(0, "5")
        entry_cell_spec_lim.delete(0, tk.END)
        entry_cell_spec_lim.insert(0, "1")
        # 隐藏 marker CCS 输入框
        label_marker_ccs.grid_remove()
        entry_marker_ccs.grid_remove()
    update_mode_params()

# ------------------------------------------------------------
# 构建图形界面
# ------------------------------------------------------------
root = tk.Tk()
status_label = ttk.Label(root, text="Ready", relief="sunken", anchor="w")
status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=2)

root.title("MPCell-Assisted Single-Cell Metabolomics")
root.geometry("1000x850")
root.minsize(500, 500)
root.resizable(True, True)

notebook = ttk.Notebook(root)
notebook.pack(fill='both', expand=True, padx=10, pady=10)

# ---------- 第一页：文件设置 ----------
frame_files = ttk.Frame(notebook)
notebook.add(frame_files, text="File Settings")

row = 0
ttk.Label(frame_files, text="MS:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
mode_var = tk.StringVar(value="no_tims")
frame_mode_sel = ttk.Frame(frame_files)
frame_mode_sel.grid(row=row, column=1, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_mode_sel, text="no_IM", variable=mode_var, value="no_tims", command=on_mode_change).pack(side='left', padx=2)
ttk.Radiobutton(frame_mode_sel, text="IM-based", variable=mode_var, value="tims", command=on_mode_change).pack(side='left', padx=2)
row += 1

ttk.Label(frame_files, text="Input file (.mzML):").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
entry_input_file = ttk.Entry(frame_files, width=80)
entry_input_file.grid(row=row, column=1, padx=(0,5), pady=5, sticky='w')
btn_browse_in = ttk.Button(frame_files, text="Browse", command=lambda: browse_input_file(entry_input_file))
btn_browse_in.grid(row=row, column=2, padx=5, pady=5)
row += 1

ttk.Label(frame_files, text="Output direction:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
entry_output_dir = ttk.Entry(frame_files, width=80)
entry_output_dir.grid(row=row, column=1, padx=(0,5), pady=5, sticky='w')
btn_browse_out = ttk.Button(frame_files, text="Browse", command=lambda: browse_output_dir(entry_output_dir))
btn_browse_out.grid(row=row, column=2, padx=5, pady=5)
row += 1

frame_files.columnconfigure(1, weight=1)

# ---------- 第二页：基本参数 ----------
frame_basic = ttk.Frame(notebook)
notebook.add(frame_basic, text="Basic Parameters")

row = 0
ttk.Label(frame_basic, text="peak_intensity_threshold:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_peak_intensity = ttk.Entry(frame_basic, width=20)
entry_peak_intensity.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_peak_intensity.insert(0, "1e4")
row += 1

label_ppm = ttk.Label(frame_basic, text="mz tolerance [ppm]:")
label_ppm.grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_ppm = ttk.Entry(frame_basic, width=20)
entry_ppm.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_ppm.insert(0, "10")
row += 1

label_ccsppm = ttk.Label(frame_basic, text="ccs tolerance:")
label_ccsppm.grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_ccsppm = ttk.Entry(frame_basic, width=20)
entry_ccsppm.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_ccsppm.insert(0, "0.05")
label_ccsppm.grid_remove()
entry_ccsppm.grid_remove()
row += 1

ttk.Label(frame_basic, text="RT unit in mzML file:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
orig_rt_unit_var = tk.StringVar(value="minute")
frame_orig_unit = ttk.Frame(frame_basic)
frame_orig_unit.grid(row=row, column=1, sticky='w', padx=0, pady=3)
ttk.Radiobutton(frame_orig_unit, text="minute", variable=orig_rt_unit_var, value="minute").pack(side='left', padx=2)
ttk.Radiobutton(frame_orig_unit, text="second", variable=orig_rt_unit_var, value="second").pack(side='left', padx=2)
row += 1

ttk.Label(frame_basic, text="rt range [min]:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
frame_rt = ttk.Frame(frame_basic)
frame_rt.grid(row=row, column=1, sticky='w', padx=0, pady=3)
entry_rt_min = ttk.Entry(frame_rt, width=8)
entry_rt_min.pack(side='left', padx=0)
entry_rt_min.insert(0, "0")
ttk.Label(frame_rt, text="~").pack(side='left', padx=2)
entry_rt_max = ttk.Entry(frame_rt, width=8)
entry_rt_max.pack(side='left', padx=0)
entry_rt_max.insert(0, "100")
row += 1

ttk.Label(frame_basic, text="Plot:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
plot_var = tk.StringVar(value="none")
frame_plot = ttk.Frame(frame_basic)
frame_plot.grid(row=row, column=1, columnspan=2, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_plot, text="None", variable=plot_var, value="none").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="Positive peak only", variable=plot_var, value="fp").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="Trash peak only", variable=plot_var, value="trash").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="All", variable=plot_var, value="all").pack(side='left', padx=2)
row += 1

ttk.Label(frame_basic, text="Plot RT range [min]:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
frame_rt_plot = ttk.Frame(frame_basic)
frame_rt_plot.grid(row=row, column=1, sticky='w', padx=0, pady=3)
entry_plot_rt_min = ttk.Entry(frame_rt_plot, width=8)
entry_plot_rt_min.pack(side='left', padx=0)
entry_plot_rt_min.insert(0, "5")
ttk.Label(frame_rt_plot, text="~").pack(side='left', padx=2)
entry_plot_rt_max = ttk.Entry(frame_rt_plot, width=8)
entry_plot_rt_max.pack(side='left', padx=0)
entry_plot_rt_max.insert(0, "10")
row += 1

ttk.Label(frame_basic, text="cell_spec:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_cell_spec = ttk.Entry(frame_basic, width=20)
entry_cell_spec.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_cell_spec.insert(0, "5")
row += 1

ttk.Label(frame_basic, text="cell_spec_lim:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_cell_spec_lim = ttk.Entry(frame_basic, width=20)
entry_cell_spec_lim.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_cell_spec_lim.insert(0, "1")
row += 1

label_eim_span = ttk.Label(frame_basic, text="eim_span:")
label_eim_span.grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_eim_span = ttk.Entry(frame_basic, width=20)
entry_eim_span.grid(row=row, column=1, padx=(0,5), pady=3, sticky='w')
entry_eim_span.insert(0, "27")
label_eim_span.grid_remove()
entry_eim_span.grid_remove()
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

ttk.Label(frame_basic, text="Filter mode:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
filter_mode_var = tk.StringVar(value="score")
frame_mode = ttk.Frame(frame_basic)
frame_mode.grid(row=row, column=1, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_mode, text="SNR-based", variable=filter_mode_var, value="snr", command=update_threshold_label).pack(side='left', padx=2)
ttk.Radiobutton(frame_mode, text="Score-based", variable=filter_mode_var, value="score", command=update_threshold_label).pack(side='left', padx=2)
row += 1

ttk.Label(frame_basic, text="Threshold:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_threshold = ttk.Entry(frame_basic, width=10)
entry_threshold.grid(row=row, column=1, sticky='w', padx=(0,2), pady=3)
threshold_hint = ttk.Label(frame_basic, text="(Score threshold, default 10)", foreground="gray")
threshold_hint.grid(row=row, column=2, sticky='w', padx=2, pady=3)
entry_threshold.insert(0, "10")
entry_threshold.bind('<KeyRelease>', on_threshold_edit)
threshold_is_default = True
update_threshold_label()

frame_basic.columnconfigure(1, weight=1)

# ---------- 第三页：离子模式标记 ----------
frame_ion = ttk.Frame(notebook)
notebook.add(frame_ion, text="Ion Mode Polarity")

# 配置列权重，确保适当拉伸
frame_ion.columnconfigure(0, weight=0, minsize=100)  # 标签列
frame_ion.columnconfigure(1, weight=1)               # 输入框列
frame_ion.columnconfigure(2, weight=0, minsize=80)   # CCS标签列
frame_ion.columnconfigure(3, weight=1)               # CCS输入框列
frame_ion.columnconfigure(4, weight=0)               # 添加行按钮列

row = 0
ttk.Label(frame_ion, text="Ion mode:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
frame_mode_ion = ttk.Frame(frame_ion)
frame_mode_ion.grid(row=row, column=1, columnspan=3, sticky='w', padx=0, pady=5)  # 扩大跨列
ion_mode_var = tk.StringVar(value="Positive")
ttk.Radiobutton(frame_mode_ion, text="Positive", variable=ion_mode_var, value="Positive", command=update_mode_params).pack(side='left', padx=5)
ttk.Radiobutton(frame_mode_ion, text="Negative", variable=ion_mode_var, value="Negative", command=update_mode_params).pack(side='left', padx=5)
row += 1

# marker 行：m/z 和 CCS 并排
ttk.Label(frame_ion, text="Marker m/z:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_marker_mz = ttk.Entry(frame_ion, width=15)
entry_marker_mz.grid(row=row, column=1, sticky='w', padx=(0,5), pady=3)

label_marker_ccs = ttk.Label(frame_ion, text="Marker CCS:")
label_marker_ccs.grid(row=row, column=2, sticky='e', padx=(0,2), pady=3)
entry_marker_ccs = ttk.Entry(frame_ion, width=15)
entry_marker_ccs.grid(row=row, column=3, sticky='w', padx=(0,5), pady=3)
row += 1

ttk.Label(frame_ion, text="intensity threshold:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=3)
entry_intensity_thres = ttk.Entry(frame_ion, width=20)
entry_intensity_thres.grid(row=row, column=1, columnspan=3, sticky='w', padx=(0,5), pady=3)
row += 1

# 背景标记区域
label_bg_hint = ttk.Label(frame_ion, text="Background markers (m/z, Intensity max):")
label_bg_hint.grid(row=row, column=0, sticky='ne', padx=(0,2), pady=5)
bg_rows_container = tk.Frame(frame_ion)
bg_rows_container.grid(row=row, column=1, columnspan=3, sticky='w', pady=5)

# 添加行按钮
btn_frame_bg = ttk.Frame(frame_ion)
btn_frame_bg.grid(row=row, column=4, rowspan=2, padx=5, sticky='n')
ttk.Button(btn_frame_bg, text="Add row", command=lambda: add_bg_row()).pack(pady=2)
row += 1

frame_ion.rowconfigure(row, weight=1)

# ---------- 第四页：运行输出 ----------
frame_output = ttk.Frame(notebook)
notebook.add(frame_output, text="Output")

result_text = tk.Text(frame_output, wrap=tk.WORD, height=5)
scrollbar = ttk.Scrollbar(frame_output, command=result_text.yview)
result_text.configure(yscrollcommand=scrollbar.set)
result_text.pack(side='left', fill='both', expand=True, padx=5, pady=5)
scrollbar.pack(side='right', fill='y')

# ---------- 第五页：强度矩阵输出 ----------
frame_quant = ttk.Frame(notebook)
notebook.add(frame_quant, text="Quantification")

row = 0
ttk.Label(frame_quant, text="Select peak_list file (.xlsx):").grid(row=row, column=0, sticky='w', padx=5, pady=5)
entry_mz_list = ttk.Entry(frame_quant, width=60)
entry_mz_list.grid(row=row, column=1, padx=(0,5), pady=5, sticky='w')
btn_browse_mz = ttk.Button(frame_quant, text="Browse", command=lambda: browse_mz_list(entry_mz_list))
btn_browse_mz.grid(row=row, column=2, padx=5, pady=5)
row += 1
ttk.Label(frame_quant, text="Note: no_IM: [ mz ]    IM-based: [ ccs , mz ]").grid(row=row, column=0, sticky='w', padx=5, pady=2)
row += 1

# 运行和终止按钮
btn_frame_quant = ttk.Frame(frame_quant)
btn_frame_quant.grid(row=row, column=0, columnspan=3, pady=10)
btn_quant_run = ttk.Button(btn_frame_quant, text="Run", command=run_quant_subprocess, width=12)
btn_quant_run.pack(side='left', padx=5)
btn_quant_terminate = ttk.Button(btn_frame_quant, text="Terminate", command=terminate_quant_subprocess, width=12, state="disabled")
btn_quant_terminate.pack(side='left', padx=5)
row += 1

# 输出文本框
ttk.Label(frame_quant, text="Quantification Output:").grid(row=row, column=0, sticky='w', padx=5, pady=2)
row += 1
quant_result_text = tk.Text(frame_quant, wrap=tk.WORD, height=20)
scrollbar_quant = ttk.Scrollbar(frame_quant, command=quant_result_text.yview)
quant_result_text.configure(yscrollcommand=scrollbar_quant.set)
quant_result_text.grid(row=row, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)
scrollbar_quant.grid(row=row, column=3, sticky='ns')
frame_quant.columnconfigure(1, weight=1)
#frame_quant.rowconfigure(row, weight=1)

# ---------- 第六页：代谢物注释 ----------
frame_annotation = ttk.Frame(notebook)
notebook.add(frame_annotation, text="Metabolite Annotation")

# 整体分为两列，用于左右并排
frame_annotation.columnconfigure(0, weight=1)
frame_annotation.columnconfigure(1, weight=1)

row = 0

# 离子模式选择（横跨两列）
ttk.Label(frame_annotation, text="Ion mode:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
frame_ion_anno = ttk.Frame(frame_annotation)
frame_ion_anno.grid(row=row, column=1, sticky='w', padx=0, pady=5, columnspan=1)
ion_mode_var_anno = tk.StringVar(value="Positive")
ttk.Radiobutton(frame_ion_anno, text="Positive", variable=ion_mode_var_anno, value="Positive").pack(side='left', padx=2)
ttk.Radiobutton(frame_ion_anno, text="Negative", variable=ion_mode_var_anno, value="Negative").pack(side='left', padx=2)
row += 1

# 待指认文件（横跨两列）
ttk.Label(frame_annotation, text="Input m/z file (.xlsx):").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
entry_mz_file = ttk.Entry(frame_annotation, width=60)
entry_mz_file.grid(row=row, column=1, padx=(0,5), pady=5, sticky='w')
btn_browse_mz = ttk.Button(frame_annotation, text="Browse", command=lambda: browse_mz_file(entry_mz_file))
btn_browse_mz.grid(row=row, column=2, padx=5, pady=5)
row += 1

# 内嵌数据库选择（横跨两列）
ttk.Label(frame_annotation, text="Embedded databases:").grid(row=row, column=0, sticky='ne', padx=(0,2), pady=5)
frame_embedded = ttk.Frame(frame_annotation)
frame_embedded.grid(row=row, column=1, columnspan=2, sticky='w', padx=(0,2), pady=5)
embedded_names = ["lcms_mda", "lcms_hela", "lcms_ninecells", "sc_paper", "hmdb", "lipidmaps"]
for i, name in enumerate(embedded_names):
    var = tk.BooleanVar(value=True)
    embedded_db_vars[name] = var
    chk = ttk.Checkbutton(frame_embedded, text=name, variable=var, command=update_db_priority_list)
    chk.grid(row=i//3, column=i%3, sticky='w', padx=5, pady=2)
row += 1

# 自定义数据库（横跨两列）
ttk.Label(frame_annotation, text="Custom database:").grid(row=row, column=0, sticky='e', padx=(0,2), pady=5)
entry_db_file = ttk.Entry(frame_annotation, width=50)
entry_db_file.grid(row=row, column=1, padx=(0,5), pady=5, sticky='w')
btn_add_custom = ttk.Button(frame_annotation, text="Add", command=add_custom_db, width=6)
btn_add_custom.grid(row=row, column=2, padx=2, pady=5, sticky='w')
btn_remove_custom = ttk.Button(frame_annotation, text="Remove", command=remove_custom_db, width=10)
btn_remove_custom.grid(row=row, column=3, padx=2, pady=5, sticky='w')
row += 1

# ========== 并列区域：左侧数据库优先级，右侧加合物优先级 ==========
priority_frame = ttk.Frame(frame_annotation)
priority_frame.grid(row=row, column=0, columnspan=4, sticky='nsew', pady=5)
priority_frame.columnconfigure(0, weight=1)
priority_frame.columnconfigure(1, weight=1)

# ----- 左侧：数据库优先级 -----
left_frame = ttk.Frame(priority_frame, relief='ridge', borderwidth=1)
left_frame.grid(row=0, column=0, sticky='nsew', padx=5)
ttk.Label(left_frame, text="Database priority", font=('Arial', 10, 'bold')).pack(pady=2)

db_container = ttk.Frame(left_frame)
db_container.pack(fill='both', expand=True, padx=5, pady=5)
listbox_db_priority = tk.Listbox(db_container, height=8, width=30, selectmode=tk.SINGLE)
listbox_db_priority.pack(side='left', fill='y')
scrollbar_db = ttk.Scrollbar(db_container, orient='vertical', command=listbox_db_priority.yview)
listbox_db_priority.configure(yscrollcommand=scrollbar_db.set)
scrollbar_db.pack(side='left', fill='y')
btn_frame_db = ttk.Frame(db_container)
btn_frame_db.pack(side='right', padx=5)
ttk.Button(btn_frame_db, text="Move up", command=move_db_up, width=10).pack(pady=2)
ttk.Button(btn_frame_db, text="Move down", command=move_db_down, width=10).pack(pady=2)
ttk.Button(btn_frame_db, text="Add", command=add_db, width=10).pack(pady=2)
ttk.Button(btn_frame_db, text="Remove", command=remove_db, width=10).pack(pady=2)

# 初始化数据库优先级列表
db_priority_list = [name for name in embedded_names]
refresh_db_listbox()

# 右侧：加合物优先级
right_frame = ttk.Frame(priority_frame, relief='ridge', borderwidth=1)
right_frame.grid(row=0, column=1, sticky='nsew', padx=5)
ttk.Label(right_frame, text="Adduct priority", font=('Arial', 10, 'bold')).pack(pady=2)
def update_adduct_listbox():
    """刷新加合物优先级列表框显示"""
    listbox_adduct_priority.delete(0, tk.END)
    for ad in adduct_priority_list:
        listbox_adduct_priority.insert(tk.END, ad)

def add_adduct():
    """从预定义加合物中选择添加"""
    ion_mode = ion_mode_var_anno.get()
    if ion_mode == "Positive":
        available = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+','[M+H-H2O]+','[M+H-2H2O]+','[M+2Na-H]+','[M+2K-H]+','[M+2Na]2+',
                     '[M+2K]2+','[M+2H]2+','[M+3H]3+','[M+4H]4+','[M.]+']
    else:
        available = ['[M-H]-','[M+Cl]-','[M+Na-2H]-','[M+OAc]-','[M+HCOO]-','[M+HCO3]-','[M-2H]2-',
                     '[M-3H]3-','[M-4H]4-','[M-CH3]-','[M.]-']
    to_add = [ad for ad in available if ad not in adduct_priority_list]
    if not to_add:
        messagebox.showinfo("Info", "All adducts have been added")
        return
    win = tk.Toplevel(root)
    win.title("Select Adduct")
    win.geometry("300x250")
    win.transient(root)
    win.grab_set()
    
    # 主框架
    main_frame = ttk.Frame(win)
    main_frame.pack(fill='both', expand=True)
    
    # 列表框
    lb = tk.Listbox(main_frame, selectmode=tk.SINGLE, height=8)
    lb.pack(fill='both', expand=True, padx=10, pady=10)
    for ad in to_add:
        lb.insert(tk.END, ad)
    
    # 按钮框架固定在底部
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(side='bottom', pady=5)
    def confirm():
        sel = lb.curselection()
        if sel:
            selected = lb.get(sel[0])
            adduct_priority_list.append(selected)
            update_adduct_listbox()
        win.destroy()
    ttk.Button(btn_frame, text="OK", command=confirm).pack()
    
def remove_adduct():
    """移除选中的加合物"""
    sel = listbox_adduct_priority.curselection()
    if sel:
        idx = sel[0]
        adduct_priority_list.pop(idx)
        update_adduct_listbox()


def move_adduct_up():
    sel = listbox_adduct_priority.curselection()
    if sel:
        idx = sel[0]
        if idx > 0:
            adduct_priority_list[idx], adduct_priority_list[idx-1] = adduct_priority_list[idx-1], adduct_priority_list[idx]
            update_adduct_listbox()
            listbox_adduct_priority.selection_set(idx-1)

def move_adduct_down():
    sel = listbox_adduct_priority.curselection()
    if sel:
        idx = sel[0]
        if idx < len(adduct_priority_list)-1:
            adduct_priority_list[idx], adduct_priority_list[idx+1] = adduct_priority_list[idx+1], adduct_priority_list[idx]
            update_adduct_listbox()
            listbox_adduct_priority.selection_set(idx+1)

def get_adduct_priority():
    return adduct_priority_list

# 初始化加合物优先级列表
def init_adduct_priority(ion_mode):
    global adduct_priority_list
    if ion_mode == "Positive":
        default = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+','[M+H-H2O]+','[M+H-2H2O]+','[M+2Na-H]+','[M+2K-H]+']
    else:
        default = ['[M-H]-','[M+Cl]-','[M+Na-2H]-','[M+OAc]-','[M+HCOO]-','[M+HCO3]-']
    adduct_priority_list = default.copy()
    update_adduct_listbox()

def on_ion_mode_change(*args):
    init_adduct_priority(ion_mode_var_anno.get())
    
adduct_container = ttk.Frame(right_frame)
adduct_container.pack(fill='both', expand=True, padx=5, pady=5)
listbox_adduct_priority = tk.Listbox(adduct_container, height=8, width=30, selectmode=tk.SINGLE)
listbox_adduct_priority.pack(side='left', fill='y')
scrollbar_adduct = ttk.Scrollbar(adduct_container, orient='vertical', command=listbox_adduct_priority.yview)
listbox_adduct_priority.configure(yscrollcommand=scrollbar_adduct.set)
scrollbar_adduct.pack(side='left', fill='y')
btn_frame_adduct = ttk.Frame(adduct_container)
btn_frame_adduct.pack(side='right', padx=5)
ttk.Button(btn_frame_adduct, text="Move up", command=move_adduct_up, width=10).pack(pady=2)
ttk.Button(btn_frame_adduct, text="Move down", command=move_adduct_down, width=10).pack(pady=2)
ttk.Button(btn_frame_adduct, text="Add", command=add_adduct, width=10).pack(pady=2)
ttk.Button(btn_frame_adduct, text="Remove", command=remove_adduct, width=10).pack(pady=2)

# 初始化加合物优先级列表
ion_mode_var_anno.trace('w', on_ion_mode_change)
init_adduct_priority(ion_mode_var_anno.get())

row += 1

# ppm 容差和最多候选数放在同一行，靠左紧凑排列
param_frame = ttk.Frame(frame_annotation)
param_frame.grid(row=row, column=0, columnspan=4, sticky='w', padx=5, pady=5)

ttk.Label(param_frame, text="m/z tolerance:").pack(side='left', padx=(0,2))
entry_ppm_anno = ttk.Entry(param_frame, width=10)
entry_ppm_anno.pack(side='left', padx=(0,10))
entry_ppm_anno.insert(0, "10")

ttk.Label(param_frame, text="Max candidates:").pack(side='left', padx=(0,2))
entry_max_candidates = ttk.Entry(param_frame, width=10)
entry_max_candidates.pack(side='left')
entry_max_candidates.insert(0, "5")

row += 1

# 运行和终止按钮
btn_frame_anno = ttk.Frame(frame_annotation)
btn_frame_anno.grid(row=row, column=0, columnspan=4, pady=10)
btn_anno_run = ttk.Button(btn_frame_anno, text="Run", command=run_annotation_subprocess, width=12)
btn_anno_run.pack(side='left', padx=5)
btn_anno_terminate = ttk.Button(btn_frame_anno, text="Terminate", command=terminate_annotation_subprocess, width=12, state="disabled")
btn_anno_terminate.pack(side='left', padx=5)
row += 1

# 输出文本框
ttk.Label(frame_annotation, text="Annotation Output:").grid(row=row, column=0, sticky='w', padx=5, pady=2)
row += 1
anno_result_text = tk.Text(frame_annotation, wrap=tk.WORD, height=10)
scrollbar_anno = ttk.Scrollbar(frame_annotation, command=anno_result_text.yview)
anno_result_text.configure(yscrollcommand=scrollbar_anno.set)
anno_result_text.grid(row=row, column=0, columnspan=4, sticky='nsew', padx=5, pady=5)
scrollbar_anno.grid(row=row, column=4, sticky='ns')
frame_annotation.columnconfigure(1, weight=1)
#frame_annotation.rowconfigure(row, weight=1)

# 初始化默认行
update_mode_params()
on_mode_change()

# ---------- 底部按钮栏 ----------
button_frame = ttk.Frame(root)
button_frame.pack(side='bottom', fill='x', pady=10)

btn_run = ttk.Button(button_frame, text="Run", command=run_subprocess_threaded, width=12)
btn_terminate = ttk.Button(button_frame, text="Terminate", command=terminate_subprocess, width=12, state="disabled")
btn_quit = ttk.Button(button_frame, text="Exit", command=root.destroy, width=12)

# 初始只打包退出按钮（固定右侧），运行/终止暂不打包
btn_quit.pack(side='right', padx=5)

def on_tab_changed(event):
    selected_tab = notebook.select()
    tab_text = notebook.tab(selected_tab, "text")
    if tab_text in ["File Settings", "Basic Parameters", "Ion Mode Polarity", "Output"]:
        # 显示运行和终止按钮（左侧）
        btn_run.pack(side='left', padx=5)
        btn_terminate.pack(side='left', padx=5)
    else:
        # 隐藏运行和终止按钮
        btn_run.pack_forget()
        btn_terminate.pack_forget()
    # 退出按钮始终存在，无需操作

notebook.bind("<<NotebookTabChanged>>", on_tab_changed)
on_tab_changed(None)   # 初始化当前页面状态

root.mainloop()


# In[ ]:




