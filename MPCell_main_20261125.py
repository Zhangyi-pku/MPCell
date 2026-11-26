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
import csv
import shutil
import tempfile

# ============================================================
# 全局状态变量
# ============================================================
current_proc = None
quant_proc = None
anno_proc = None
metcell_proc = None
metcell_container_name = None
metcell_volume_name = None
worker_tims_proc = None
btn_worker_tims_run = None
# ---- no_tims 模式：Parameter Setting 页使用的全局变量 ----
bg_rows = []                 # 背景标记行 [{frame, mz, intensity}]
bg_rows_container = None     # 背景标记行容器
entry_peak_intensity = None
entry_ppm = None
entry_rt_min = None
entry_rt_max = None
entry_marker_mz = None
entry_intensity_thres = None
entry_cell_eic_span = None
entry_spec_num_lim = None
entry_cell_eic_span_lim = None
entry_cell_occurrence = None
entry_threshold = None
threshold_hint = None
threshold_is_default = True
plot_var = None
entry_plot_rt_min = None
entry_plot_rt_max = None
ion_mode_var = None
# ---- Cell-specific Peak Parameters 页变量 ----
var_use_snr = None
var_use_match = None
var_use_rsd = None
entry_snr_thres = None
entry_snr_lower = None
entry_match_thres = None
entry_match_lower = None
entry_rsd_thres = None
entry_rsd_lower = None
cell_specific_sf = None

# ---- tims 模式：MetCell EIC 页全部变量放在 mc_vars ----
mc_vars = {}
mc_bg_rows = []              # MetCell 背景标记行 [{frame, mz, ccs, intensity}]
mc_bg_rows_container = None
time_seg_rows = []
seg_rows_container = None

# ---- Metabolite Annotation 状态 ----
embedded_db_vars = {}
db_priority_list = []
custom_db_path = None
adduct_priority_list = []

# ---- 各页 Frame 引用 ----
file_settings_sf = None
param_setting_sf = None
output_sf = None
quant_sf = None
anno_sf = None
metcell_sf = None
mc_output_sf = None
cell_specific_tims_sf = None
cell_specific_tims_sf = None
btn_cs_tims_terminate = None   # 新页 Terminate

# ---- 镜像候选 ----
METCELL_IMAGE_CANDIDATES = [
    "zhulab/metcell-r:latest",
    "docker.1ms.run/zhulab/metcell-r:latest",
    "docker.m.daocloud.io/zhulab/metcell-r:latest",
]


# ============================================================
# 可滚动容器
# ============================================================
class ScrollableFrame(ttk.Frame):
    """一个带纵向滚动条的 ttk.Frame 包装。"""
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        try:
            bg = ttk.Style().lookup('TFrame', 'background') or '#f0f0f0'
        except Exception:
            bg = '#f0f0f0'

        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0, background=bg)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.inner = ttk.Frame(self.canvas)
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")

        self.inner.bind("<Configure>", self._on_inner_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)
        self.canvas.bind("<Enter>", self._bind_mousewheel)
        self.canvas.bind("<Leave>", self._unbind_mousewheel)

    def _on_inner_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self.canvas.itemconfigure(self._win, width=event.width)

    def _bind_mousewheel(self, event):
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.canvas.bind_all("<Button-4>", self._on_mousewheel)
        self.canvas.bind_all("<Button-5>", self._on_mousewheel)

    def _unbind_mousewheel(self, event):
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")

    def _on_mousewheel(self, event):
        w = event.widget
        if isinstance(w, (tk.Text, tk.Listbox, tk.Canvas, ttk.Treeview)):
            return
        if event.num == 4:
            units = -1
        elif event.num == 5:
            units = 1
        else:
            delta = event.delta
            if delta == 0:
                return
            if abs(delta) >= 120:
                units = -int(delta / 120)
            else:
                units = -1 if delta > 0 else 1
        self.canvas.yview_scroll(units, "units")


# ============================================================
# 文件/文件夹浏览回调
# ============================================================
def browse_input_file(entry_widget):
    path = filedialog.askopenfilename(
        title="Select Input File",
        filetypes=[
            ("MS files", "*.mzxml *.mzML"),
            ("mzXML files", "*.mzxml"),
            ("mzML files", "*.mzml"),
            ("All files", "*.*"),
        ]
    )
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)
        if not entry_output_dir.get().strip():
            entry_output_dir.delete(0, tk.END)
            entry_output_dir.insert(0, get_default_output_dir(path))


def browse_output_dir(entry_widget):
    path = filedialog.askdirectory(title="Select Output Folder")
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)
        if not entry_mz_list.get().strip():
            entry_mz_list.delete(0, tk.END)
            entry_mz_list.insert(0, os.path.join(path, "grouping", "RT_del", "c_peak_fp.xlsx"))


def browse_mz_list(entry_widget):
    path = filedialog.askopenfilename(
        title="Select peak_list Excel File",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
    )
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)


def browse_mz_file(entry_widget):
    path = filedialog.askopenfilename(
        title="Select m/z Input File (Excel)",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
    )
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)


def get_default_output_dir(input_file):
    if not input_file:
        return ""
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    return os.path.join(os.path.dirname(input_file), base_name)


# ============================================================
# 数据库优先级（no_tims / Metabolite Annotation）
# ============================================================
def update_db_priority_list():
    global db_priority_list, custom_db_path
    new_embedded = [name for name, var in embedded_db_vars.items() if var.get()]
    old_order = db_priority_list.copy()
    result = []
    for db in old_order:
        if db == "Custom":
            if custom_db_path:
                result.append(db)
        elif db in new_embedded:
            result.append(db)
    for db in new_embedded:
        if db not in result:
            result.append(db)
    if custom_db_path and "Custom" not in result:
        result.append("Custom")
    db_priority_list = result
    refresh_db_listbox()


def refresh_db_listbox():
    listbox_db_priority.delete(0, tk.END)
    for db in db_priority_list:
        listbox_db_priority.insert(tk.END, db)


def add_custom_db():
    global custom_db_path
    path = filedialog.askopenfilename(
        title="Select Custom Database",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
    )
    if path:
        custom_db_path = path
        entry_db_file.delete(0, tk.END)
        entry_db_file.insert(0, path)
        update_db_priority_list()


def remove_custom_db():
    global custom_db_path
    custom_db_path = None
    entry_db_file.delete(0, tk.END)
    update_db_priority_list()
    listbox_db_priority.selection_clear(0, tk.END)


def add_db():
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
    sel = listbox_db_priority.curselection()
    if not sel:
        messagebox.showinfo("Info", "Please select a database to remove")
        return
    idx = sel[0]
    removed = db_priority_list.pop(idx)
    if removed in embedded_db_vars:
        embedded_db_vars[removed].set(False)
    if removed == "Custom":
        global custom_db_path
        custom_db_path = None
        entry_db_file.delete(0, tk.END)
    refresh_db_listbox()
    listbox_db_priority.selection_clear(0, tk.END)


def move_db_up():
    sel = listbox_db_priority.curselection()
    if not sel:
        return
    idx = sel[0]
    if idx > 0:
        db_priority_list[idx], db_priority_list[idx - 1] = db_priority_list[idx - 1], db_priority_list[idx]
        refresh_db_listbox()
        listbox_db_priority.selection_set(idx - 1)


def move_db_down():
    sel = listbox_db_priority.curselection()
    if not sel:
        return
    idx = sel[0]
    if idx < len(db_priority_list) - 1:
        db_priority_list[idx], db_priority_list[idx + 1] = db_priority_list[idx + 1], db_priority_list[idx]
        refresh_db_listbox()
        listbox_db_priority.selection_set(idx + 1)


# ============================================================
# 加合物优先级（no_tims / Metabolite Annotation）
# ============================================================
def update_adduct_listbox():
    listbox_adduct_priority.delete(0, tk.END)
    for ad in adduct_priority_list:
        listbox_adduct_priority.insert(tk.END, ad)


def add_adduct():
    ion_mode = ion_mode_var_anno.get()
    if ion_mode == "Positive":
        available = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+', '[M+H-H2O]+', '[M+H-2H2O]+',
                     '[M+2Na-H]+', '[M+2K-H]+', '[M+2Na]2+', '[M+2K]2+', '[M+2H]2+',
                     '[M+3H]3+', '[M+4H]4+', '[M.]+']
    else:
        available = ['[M-H]-', '[M+Cl]-', '[M+Na-2H]-', '[M+OAc]-', '[M+HCOO]-', '[M+HCO3]-',
                     '[M-2H]2-', '[M-3H]3-', '[M-4H]4-', '[M-CH3]-', '[M.]-']
    to_add = [ad for ad in available if ad not in adduct_priority_list]
    if not to_add:
        messagebox.showinfo("Info", "All adducts have been added")
        return
    win = tk.Toplevel(root)
    win.title("Select Adduct")
    win.geometry("300x250")
    win.transient(root)
    win.grab_set()
    main_frame = ttk.Frame(win)
    main_frame.pack(fill='both', expand=True)
    lb = tk.Listbox(main_frame, selectmode=tk.SINGLE, height=8)
    lb.pack(fill='both', expand=True, padx=10, pady=10)
    for ad in to_add:
        lb.insert(tk.END, ad)
    btn_frame = ttk.Frame(main_frame)
    btn_frame.pack(side='bottom', pady=5)

    def confirm():
        sel = lb.curselection()
        if sel:
            adduct_priority_list.append(lb.get(sel[0]))
            update_adduct_listbox()
        win.destroy()
    ttk.Button(btn_frame, text="OK", command=confirm).pack()


def remove_adduct():
    sel = listbox_adduct_priority.curselection()
    if sel:
        adduct_priority_list.pop(sel[0])
        update_adduct_listbox()


def move_adduct_up():
    sel = listbox_adduct_priority.curselection()
    if sel and sel[0] > 0:
        i = sel[0]
        adduct_priority_list[i], adduct_priority_list[i - 1] = adduct_priority_list[i - 1], adduct_priority_list[i]
        update_adduct_listbox()
        listbox_adduct_priority.selection_set(i - 1)


def move_adduct_down():
    sel = listbox_adduct_priority.curselection()
    if sel and sel[0] < len(adduct_priority_list) - 1:
        i = sel[0]
        adduct_priority_list[i], adduct_priority_list[i + 1] = adduct_priority_list[i + 1], adduct_priority_list[i]
        update_adduct_listbox()
        listbox_adduct_priority.selection_set(i + 1)


def get_adduct_priority():
    return adduct_priority_list


def init_adduct_priority(ion_mode):
    global adduct_priority_list
    if ion_mode == "Positive":
        default = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+', '[M+H-H2O]+', '[M+H-2H2O]+', '[M+2Na-H]+', '[M+2K-H]+']
    else:
        default = ['[M-H]-', '[M+Cl]-', '[M+Na-2H]-', '[M+OAc]-', '[M+HCOO]-', '[M+HCO3]-']
    adduct_priority_list = default.copy()
    update_adduct_listbox()


def on_ion_mode_change(*args):
    init_adduct_priority(ion_mode_var_anno.get())


# ============================================================
# 命令获取
# ============================================================
def get_worker_command():
    if getattr(sys, 'frozen', False):
        return [os.path.join(os.path.dirname(sys.executable), "worker.exe")]
    return [sys.executable, "worker_20261125.py"]


def get_quant_command():
    if getattr(sys, 'frozen', False):
        return [os.path.join(os.path.dirname(sys.executable), 'quantification.exe')]
    return [sys.executable, "quantification_20261122.py"]


def get_annotation_command():
    if getattr(sys, 'frozen', False):
        return [os.path.join(os.path.dirname(sys.executable), 'annotation_worker.exe')]
    return [sys.executable, "annotation_worker_20261122.py"]

def get_worker_tims_command():
    if getattr(sys, 'frozen', False):
        base_dir = os.path.dirname(sys.executable)
        return [os.path.join(base_dir, "worker_tims.exe")]
    else:
        return [sys.executable, "worker_tims_20261125.py"]

# ============================================================
# no_tims 背景标记行
# ============================================================
def add_bg_row(mz_val="", int_val=""):
    row_frame = tk.Frame(bg_rows_container)
    row_frame.pack(fill='x', pady=2)

    row_frame.columnconfigure(2, weight=1)  # 弹簧列

    entry_mz = tk.Entry(row_frame, width=15)
    entry_mz.grid(row=0, column=0, padx=(0, 2), sticky='w')
    entry_mz.insert(0, mz_val)

    entry_int = tk.Entry(row_frame, width=12)
    entry_int.grid(row=0, column=1, padx=(0, 2), sticky='w')
    entry_int.insert(0, int_val)

    del_btn = tk.Button(row_frame, text="X", width=3,
                        bg="#ffdddd", activebackground="#ffaaaa",
                        command=lambda: delete_bg_row(row_frame))
    del_btn.grid(row=0, column=3, padx=5, sticky='e')

    bg_rows.append({'frame': row_frame, 'mz': entry_mz, 'intensity': entry_int})


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
    """no_tims 模式专用：defaults 为 [(mz, intensity), ...]"""
    clear_bg_rows()
    for mz, ints in defaults:
        add_bg_row(mz, ints)


# ============================================================
# MetCell EIC 背景标记行（独立容器 + 独立列表）
# ============================================================
def add_mc_bg_row(mz_val="", ccs_val="", int_val=""):
    row_frame = tk.Frame(mc_bg_rows_container)
    row_frame.pack(fill='x', pady=1)

    row_frame.columnconfigure(3, weight=1)

    entry_mz = tk.Entry(row_frame, width=10)
    entry_mz.grid(row=0, column=0, padx=(0, 2), sticky='w')
    entry_mz.insert(0, mz_val)

    entry_ccs = tk.Entry(row_frame, width=10)
    entry_ccs.grid(row=0, column=1, padx=(0, 2), sticky='w')
    entry_ccs.insert(0, ccs_val)

    entry_int = tk.Entry(row_frame, width=10)
    entry_int.grid(row=0, column=2, padx=(0, 2), sticky='w')
    entry_int.insert(0, int_val)

    def del_row():
        row_frame.destroy()
        for i, r in enumerate(mc_bg_rows):
            if r['frame'] is row_frame:
                mc_bg_rows.pop(i)
                break

    del_btn = tk.Button(row_frame, text="X", width=3,
                        bg="#ffdddd", activebackground="#ffaaaa",
                        command=del_row)
    del_btn.grid(row=0, column=4, padx=5, sticky='e')

    mc_bg_rows.append({'frame': row_frame, 'mz': entry_mz,
                       'ccs': entry_ccs, 'intensity': entry_int})


def clear_mc_bg_rows():
    for row in mc_bg_rows:
        row['frame'].destroy()
    mc_bg_rows.clear()

def init_mc_bg_defaults():
    clear_mc_bg_rows()
    ion_mode = mc_vars["ion_mode"].get() if "ion_mode" in mc_vars else "positive"
    if ion_mode == "positive":
        defaults = [
            ("760.5845", "1.402", "0"),
        ]
    else:
        defaults = [
            ("306.0767", "0",  "0"),
            ("804.5732", "0",  "0"),
        ]
    for mz, ccs, ints in defaults:
        add_mc_bg_row(mz, ccs, ints)


# ============================================================
# no_tims 模式：参数默认值
# ============================================================
def update_mode_params():
    """no_tims 模式下根据正负离子更新 marker m/z、强度阈值、背景行默认值"""
    if ion_mode_var is None or entry_marker_mz is None:
        return
    if ion_mode_var.get() == "Positive":
        marker_mz = "760.5851"
        intensity_thres_val = "1e6"
        bg_defaults = [("760.5851", "0"), ("732.5532", "0")]
    else:
        marker_mz = "306.0767"
        intensity_thres_val = "5e6"
        bg_defaults = [("306.0767", "2e6"), ("804.5732", "0"), ("742.5393", "0")]

    entry_marker_mz.delete(0, tk.END)
    entry_marker_mz.insert(0, marker_mz)
    entry_intensity_thres.delete(0, tk.END)
    entry_intensity_thres.insert(0, intensity_thres_val)
    populate_bg_rows(bg_defaults)


def on_threshold_edit(event):
    global threshold_is_default
    threshold_is_default = False


# ============================================================
# no_tims 参数收集（只用于 no_tims）
# ============================================================
def get_parameters():
    params = {}
    params['peak_intensity_thres'] = float(entry_peak_intensity.get())
    params['ppm'] = (1e-6) * float(entry_ppm.get())
    params['rt_range_min'] = float(entry_rt_min.get()) * 60.0
    params['rt_range_max'] = float(entry_rt_max.get()) * 60.0
    params['cell_eic_span'] = int(entry_cell_eic_span.get())
    params['cell_eic_span_lim'] = int(entry_cell_eic_span_lim.get())
    params['spec_num_lim'] = int(entry_spec_num_lim.get())
    params['data_type'] = 'no_tims'

    try:
        params['cell_occurrence_threshold'] = float(entry_cell_occurrence.get()) / 100.0
    except Exception:
        params['cell_occurrence_threshold'] = 0.05

    params['marker'] = float(entry_marker_mz.get())
    params['intensity_thres'] = float(entry_intensity_thres.get())

    bg_marker_list = []
    bg_intensity_list = []
    for row in bg_rows:
        mz_str = row['mz'].get().strip()
        int_str = row['intensity'].get().strip()
        if mz_str == "":
            continue
        bg_marker_list.append(float(mz_str))
        bg_intensity_list.append(float(int_str) if int_str else 0.0)
    params['bg_marker'] = bg_marker_list
    params['bg_intensity_max'] = bg_intensity_list

    params['input_file'] = entry_input_file.get().strip()
    output_dir = entry_output_dir.get().strip()
    if not output_dir and params['input_file']:
        output_dir = get_default_output_dir(params['input_file'])
    params['output_dir'] = output_dir
    params['plot_option'] = plot_var.get()
    params['plot_rt_min'] = float(entry_plot_rt_min.get()) * 60.0
    params['plot_rt_max'] = float(entry_plot_rt_max.get()) * 60.0
    params['filter_threshold'] = float(entry_threshold.get())
    params['use_snr']   = bool(var_use_snr.get())
    params['snr_thres'] = float(entry_snr_thres.get())
    params['snr_lower'] = float(entry_snr_lower.get())

    params['use_match']   = bool(var_use_match.get())
    params['match_thres'] = float(entry_match_thres.get())
    params['match_lower'] = float(entry_match_lower.get())

    params['use_rsd']   = bool(var_use_rsd.get())
    params['rsd_thres'] = float(entry_rsd_thres.get())
    params['rsd_lower'] = float(entry_rsd_lower.get())
    return params


# ============================================================
# 终止子进程
# ============================================================
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


# ============================================================
# no_tims 主流程子进程
# ============================================================
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
            proc = subprocess.Popen(
                get_worker_command() + [param_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, bufsize=1, env=env
            )
            current_proc = proc
            t_out = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_err = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_out.daemon = t_err.daemon = True
            t_out.start(); t_err.start()
            rc = proc.wait()
            t_out.join(timeout=1); t_err.join(timeout=1)
            output_queue.put(("DONE", rc))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    rc = item[1]
                    if rc == 0:
                        mc_result_text.insert(tk.END,
                            "\n*** worker_tims completed successfully ***\n")
                    else:
                        mc_result_text.insert(tk.END,
                            f"\n*** worker_tims FAILED (rc={rc}) ***\n")
                    try:
                        os.remove(param_file_path)
                    except Exception:
                        pass
                    globals()['worker_tims_proc'] = None
                    try:
                        if btn_worker_tims_run is not None:
                            btn_worker_tims_run.config(state="normal")
                        btn_mc_run.config(state="normal")
                        if metcell_proc is None:
                            btn_mc_terminate.config(state="disabled")
                    except Exception:
                        pass
                    mc_result_text.see(tk.END)
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


# ============================================================
# 定量子进程
# ============================================================
def run_quant_subprocess():
    global quant_proc
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

    params = {
        'mz_list': mz_list_path,
        'input_file': input_file,
        'output_dir': output_dir,
        'data_type': 'no_tims',
    }
    try:
        params['ppm'] = float(entry_ppm.get()) * 1e-6
    except Exception:
        pass

    param_file_path = os.path.join(output_dir, 'temp_quant_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Unable to create parameter file: {e}")
        return

    btn_quant_run.config(state="disabled")
    btn_quant_terminate.config(state="normal")
    quant_result_text.delete(1.0, tk.END)
    quant_result_text.insert(tk.END, "Starting quantification...\n")
    root.update()

    output_queue = queue.Queue()

    def read_stream(stream, name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global quant_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            proc = subprocess.Popen(
                get_quant_command() + [param_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, bufsize=1, env=env
            )
            quant_proc = proc
            t_out = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_err = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_out.daemon = t_err.daemon = True
            t_out.start(); t_err.start()
            rc = proc.wait()
            t_out.join(timeout=1); t_err.join(timeout=1)
            output_queue.put(("DONE", rc))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    rc = item[1]
                    if rc == 0:
                        quant_result_text.insert(tk.END, "\n*** Quantification completed successfully ***\n")
                    else:
                        quant_result_text.insert(tk.END, f"\n*** Quantification exited with error code {rc} ***\n")
                    try:
                        os.remove(param_file_path)
                    except Exception:
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


# ============================================================
# 注释子进程
# ============================================================
def run_annotation_subprocess():
    global anno_proc, db_priority_list, custom_db_path
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
    except Exception:
        ppm_tol = 10
    try:
        max_candidates = int(entry_max_candidates.get())
    except Exception:
        max_candidates = 5
    adduct_priority = get_adduct_priority()
    if not adduct_priority:
        messagebox.showwarning("Warning", "Please select at least one adduct type")
        return
    ion_mode = ion_mode_var_anno.get()
    if not db_priority_list:
        messagebox.showwarning("Warning", "Please select at least one embedded database")
        return
    custom_db = entry_db_file.get().strip()
    if custom_db and not os.path.exists(custom_db):
        messagebox.showwarning("Warning", f"Custom database file does not exist: {custom_db}")
        return

    params = {
        'input_mz_file': mz_file,
        'output_dir': output_dir,
        'ppm_tol': ppm_tol,
        'max_candidates': max_candidates,
        'adduct_priority': adduct_priority,
        'ion_mode': ion_mode,
        'db_priority': db_priority_list,
        'custom_db': custom_db,
    }
    param_file_path = os.path.join(output_dir, 'temp_annotation_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        messagebox.showerror("Parameter Error", f"Unable to create parameter file: {e}")
        return

    btn_anno_run.config(state="disabled")
    btn_anno_terminate.config(state="normal")
    anno_result_text.delete(1.0, tk.END)
    anno_result_text.insert(tk.END, "Starting annotation...\n")
    root.update()

    output_queue = queue.Queue()

    def read_stream(stream, name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global anno_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            proc = subprocess.Popen(
                get_annotation_command() + [param_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, bufsize=1, env=env
            )
            anno_proc = proc
            t_out = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_err = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_out.daemon = t_err.daemon = True
            t_out.start(); t_err.start()
            rc = proc.wait()
            t_out.join(timeout=1); t_err.join(timeout=1)
            output_queue.put(("DONE", rc))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    rc = item[1]
                    if rc == 0:
                        anno_result_text.insert(tk.END, "\n*** Annotation completed successfully ***\n")
                    else:
                        anno_result_text.insert(tk.END, f"\n*** Annotation exited with error code {rc} ***\n")
                    try:
                        os.remove(param_file_path)
                    except Exception:
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


# ============================================================
# MetCell EIC 辅助
# ============================================================
def _browse_into(var, is_dir=True):
    path = filedialog.askdirectory(title="Select Folder") if is_dir else \
           filedialog.askopenfilename(title="Select File")
    if path:
        var.set(path)


def _get_metcell_files_dir():
    if getattr(sys, 'frozen', False):
        base = os.path.dirname(sys.executable)
    else:
        base = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, 'metcell_files')


def _check_docker_available():
    try:
        r = subprocess.run(["docker", "--version"], capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=10)
        return r.returncode == 0
    except Exception:
        return False


def _list_docker_images():
    try:
        r = subprocess.run(["docker", "images", "--format", "{{.Repository}}:{{.Tag}}"],
                           capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=15)
        if r.returncode != 0:
            return []
        return [line.strip() for line in r.stdout.splitlines() if line.strip()]
    except Exception:
        return []


def _find_existing_image(candidates):
    existing = _list_docker_images()
    for name in candidates:
        if name in existing:
            return name
        last_part = name.rsplit("/", 1)[-1]
        if ":" not in last_part and (name + ":latest") in existing:
            return name + ":latest"
    return None


def _docker_run_streaming(cmd, output_queue, proc_holder):
    env = os.environ.copy()
    env['PYTHONUNBUFFERED'] = '1'
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding='utf-8', errors='replace',
                            bufsize=1, env=env)
    proc_holder['proc'] = proc

    def read_stream(stream, name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((name, f"[stream {name} read error] {e}"))

    t_out = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
    t_err = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
    t_out.daemon = t_err.daemon = True
    t_out.start(); t_err.start()
    rc = proc.wait()
    t_out.join(timeout=1); t_err.join(timeout=1)
    proc_holder['proc'] = None
    return rc


def _cleanup_volume(vol_name, container_name=None):
    if container_name:
        try:
            subprocess.run(["docker", "kill", container_name], capture_output=True,
                           encoding='utf-8', errors='replace', timeout=15)
        except Exception:
            pass
        try:
            subprocess.run(["docker", "rm", "-f", container_name], capture_output=True,
                           encoding='utf-8', errors='replace', timeout=15)
        except Exception:
            pass
    if vol_name:
        try:
            subprocess.run(["docker", "volume", "rm", "-f", vol_name], capture_output=True,
                           encoding='utf-8', errors='replace', timeout=15)
        except Exception:
            pass


# ============================================================
# MetCell 时间段行
# ============================================================
def add_time_seg_row(start="", end=""):
    row_frame = tk.Frame(seg_rows_container)
    row_frame.pack(fill='x', pady=1)

    row_frame.columnconfigure(2, weight=1)

    entry_start = tk.Entry(row_frame, width=15)
    entry_start.grid(row=0, column=0, padx=(0, 2), sticky='w')
    entry_start.insert(0, start)

    entry_end = tk.Entry(row_frame, width=15)
    entry_end.grid(row=0, column=1, padx=(0, 2), sticky='w')
    entry_end.insert(0, end)

    del_btn = tk.Button(row_frame, text="X", width=3,
                        bg="#ffdddd", activebackground="#ffaaaa",
                        command=lambda: delete_time_seg_row(row_frame))
    del_btn.grid(row=0, column=3, padx=5, sticky='e')

    time_seg_rows.append({'frame': row_frame,
                          'start': entry_start, 'end': entry_end})


def delete_time_seg_row(row_frame):
    for i, row in enumerate(time_seg_rows):
        if row['frame'] == row_frame:
            row_frame.destroy()
            time_seg_rows.pop(i)
            break


def clear_time_seg_rows():
    for row in time_seg_rows:
        row['frame'].destroy()
    time_seg_rows.clear()


def _default_mc_output_dir(d_folder):
    """根据 .d 文件夹自动推导默认输出目录：
       <.d 同级目录>/<去掉 .d 后缀的目录名>_metcell_out
    """
    if not d_folder:
        return ""
    d_folder = d_folder.rstrip('/\\')
    parent = os.path.dirname(d_folder)
    base = os.path.basename(d_folder)
    if base.lower().endswith('.d'):
        base = base[:-2]
    if not base or not parent:
        return ""
    return os.path.join(parent, base + '_metcell_out')


def update_data_file_display():
    if "d_folder" not in mc_vars:
        return
    d_folder = mc_vars["d_folder"].get().strip()
    if d_folder:
        mc_vars["data_file"].set(os.path.basename(d_folder.rstrip('/\\')))
        # 关键：只在 output_dir 为空时自动填充，不覆盖用户手改
        if "output_dir" in mc_vars and not mc_vars["output_dir"].get().strip():
            mc_vars["output_dir"].set(_default_mc_output_dir(d_folder))
    else:
        mc_vars["data_file"].set("")


# ============================================================
# MetCell 主流程
# ============================================================
def run_metcell_subprocess():
    global metcell_proc, metcell_container_name, metcell_volume_name

    d_folder = mc_vars["d_folder"].get().strip()
    out_dir = mc_vars["output_dir"].get().strip()
    image = mc_vars["docker_image"].get().strip()

    if not d_folder or not os.path.isdir(d_folder):
        messagebox.showwarning("Warning", "Please select a valid Bruker .d folder")
        return
    if not d_folder.lower().endswith('.d'):
        if not messagebox.askyesno("Confirm",
                f"The selected folder does not end with '.d':\n{d_folder}\n\nContinue anyway?"):
            return
    if not out_dir:
        messagebox.showwarning("Warning", "Please select an output directory")
        return
    if not image:
        messagebox.showwarning("Warning", "Please enter a Docker image name")
        return
    if not time_seg_rows:
        messagebox.showwarning("Warning",
            "Please add at least one time segment before running MetCell.")
        return

    segments = []
    for i, row in enumerate(time_seg_rows):
        s_text = row['start'].get().strip()
        e_text = row['end'].get().strip()
        if not s_text or not e_text:
            messagebox.showerror("Parameter error", f"Time segment #{i+1}: start_time or end_time is empty.")
            return
        try:
            s_val, e_val = float(s_text), float(e_text)
        except ValueError:
            messagebox.showerror("Parameter error", f"Time segment #{i+1}: start_time and end_time must be numeric.")
            return
        if s_val < 0 or e_val < 0:
            messagebox.showerror("Parameter error", f"Time segment #{i+1}: time cannot be negative.")
            return
        if s_val >= e_val:
            messagebox.showerror("Parameter error", f"Time segment #{i+1}: start_time must be less than end_time.")
            return
        segments.append((s_val, e_val))

    metcell_files = _get_metcell_files_dir()
    required = ["run.R", "extract_eic_full.R", "libtimsdata.so"]
    missing = [f for f in required if not os.path.exists(os.path.join(metcell_files, f))]
    if missing:
        messagebox.showerror("Missing files",
            f"The following files are missing from:\n{metcell_files}\n\n"
            + "\n".join(f"  - {f}" for f in missing))
        return

    if not _check_docker_available():
        messagebox.showerror("Docker not available",
            "Cannot run 'docker'.\n\nPlease install Docker Desktop and make sure it is running.")
        return

    candidates = []
    if image:
        candidates.append(image)
    for c in METCELL_IMAGE_CANDIDATES:
        if c not in candidates:
            candidates.append(c)

    resolved_image = _find_existing_image(candidates)

    if resolved_image is None:
        existing_metcell = [x for x in _list_docker_images() if "metcell" in x.lower()]
        detail = "\n".join(existing_metcell) if existing_metcell else "(none)"
        ans = messagebox.askyesno("Image not found",
            f"None of the following images are available locally:\n"
            + "\n".join(f"  - {c}" for c in candidates)
            + f"\n\nImages containing 'metcell' on this machine:\n{detail}\n\n"
            + "Try pulling now? (will try each name in order)")
        if not ans:
            return

        mc_result_text.delete(1.0, tk.END)
        mc_result_text.insert(tk.END, "Pulling image...\n")
        root.update()
        pulled_ok = False
        for c in candidates:
            mc_result_text.insert(tk.END, f"\n>>> docker pull {c}\n")
            root.update()
            try:
                r = subprocess.run(["docker", "pull", c], capture_output=True, text=True,
                                   encoding='utf-8', errors='replace', timeout=1800)
                if r.stdout:
                    mc_result_text.insert(tk.END, r.stdout)
                if r.stderr:
                    mc_result_text.insert(tk.END, r.stderr)
                root.update()
                if r.returncode == 0:
                    resolved_image = c
                    pulled_ok = True
                    mc_result_text.insert(tk.END, f"\n[OK] pulled {c}\n")
                    break
                else:
                    mc_result_text.insert(tk.END, f"\n[FAIL] {c}\n")
            except subprocess.TimeoutExpired:
                mc_result_text.insert(tk.END, f"\n[FAIL] {c} (timeout)\n")
            except Exception as e:
                mc_result_text.insert(tk.END, f"\n[FAIL] {c}: {e}\n")
        if not pulled_ok:
            messagebox.showerror("Pull failed",
                "All image names failed to pull.\nCheck your network / Docker login status, or pull manually.")
            return

    image = resolved_image

    d_folder_abs = os.path.abspath(d_folder)
    d_parent_abs = os.path.dirname(d_folder_abs)
    d_basename = os.path.basename(d_folder_abs)
    out_abs = os.path.abspath(out_dir)
    metcell_abs = os.path.abspath(metcell_files)
    os.makedirs(out_abs, exist_ok=True)

    tmp_csv_dir = tempfile.mkdtemp(prefix="mpcell_metcell_")
    tmp_csv_path = os.path.join(tmp_csv_dir, "time_limit_table.csv")
    try:
        with open(tmp_csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(["data_file", "start_time", "end_time"])
            for s_val, e_val in segments:
                s_out = int(s_val) if float(s_val).is_integer() else s_val
                e_out = int(e_val) if float(e_val).is_integer() else e_val
                writer.writerow([d_basename, s_out, e_out])
    except Exception as e:
        shutil.rmtree(tmp_csv_dir, ignore_errors=True)
        messagebox.showerror("Error", f"Failed to write time_limit_table.csv: {e}")
        return

    import time as _time
    ts = int(_time.time())
    vol_name = f"mpcell_metcell_vol_{ts}"
    container_name = f"mpcell_metcell_run_{ts}"
    metcell_volume_name = vol_name
    metcell_container_name = container_name

    # 注意：修正 "marker_span" -> "marker_eic_span"（避免 KeyError）
    param_map = {
        "METCELL_ION_MODE":               "ion_mode",
        "METCELL_N_SLAVES":               "n_slaves",
        "METCELL_MARKER_MZ":              "marker_mz",
        "METCELL_MARKER_MOBILITY":        "marker_mobility",
        "METCELL_MZ_TOL_PPM":             "mz_tol_ppm",
        "METCELL_MOBILITY_RANGE":         "mobility_range",
        "METCELL_INT_THRES_UPPER":        "int_upper",
        "METCELL_INT_THRES_LOWER":        "int_lower",
        "METCELL_MARKER_EIC_PEAK_SPAN":   "cell_eic_span",
        "METCELL_MARKER_EIM_PEAK_SPAN":   "marker_eim_span",
        "METCELL_QUERY_POOL_SIZE":        "query_pool",
        "METCELL_SUPERPOSITION_MZ_TOL":   "super_mz_tol",
        "METCELL_CELL_NUMBER":            "cell_number",
        "METCELL_SUPERPOSITION_POOL":     "super_pool",
        "METCELL_SEARCH_PEAK_MZ_TOL":     "search_mz_tol",
        "METCELL_SEARCH_N_SKIP":          "search_n_skip",
        "METCELL_PEAK_TARGET_LENGTH":     "peak_len",
        "METCELL_SMOOTH_WINDOW":          "smooth_window",
        "METCELL_PEAK_SPAN_DETECTION":    "span_det",
        "METCELL_PEAK_SPAN_INTEGRATION":  "span_intg",
        "METCELL_SIGNAL_SD_THRESHOLD":    "signal_sd",
        "METCELL_SINGLE_CHARGE_SLOPE":    "chg_slope",
        "METCELL_SINGLE_CHARGE_INTERCEPT": "chg_intercept",
        "METCELL_EXTRACT_MZ_TOL":         "extract_mz",
        "METCELL_MOBILITY_TOL":           "extract_mob",
        "METCELL_MIN_CELL_NUMBER":        "min_cell_num",
        "METCELL_EXTRACT_POOL":           "extract_pool",
        "METCELL_ISO_MZ_TOL":             "iso_mz_tol",
        "METCELL_ISO_MOBILITY_TOL":       "iso_mob_tol",
        "METCELL_ISO_DELTA":              "iso_delta",
        "METCELL_ISO_MAX_NUM":            "iso_max_num",
        "METCELL_ISO_INT_RATIO_CHECK":    "iso_ratio_chk",
        "METCELL_ISO_INT_CUTOFF":         "iso_ratio_cut",
        "METCELL_SC_TYPE_CCS":            "sc_type_ccs",
        "METCELL_SC_TOL_CCS":             "sc_tol_ccs",
        "METCELL_SC_TOL_MZ":              "sc_tol_mz",
        "METCELL_SC_ADDUCT":              "sc_adduct",
        "METCELL_EIC_MZ_TOL":             "eic_mz_tol",
        "METCELL_EIC_CHUNK_SIZE":         "eic_chunk",
        "METCELL_EIC_N_THREADS":          "eic_n_threads",
    }
    env_pairs = []
    for env_key, var_key in param_map.items():
        val = mc_vars[var_key].get().strip()
        if val == "":
            messagebox.showerror("Parameter error", f"Parameter '{var_key}' is empty")
            shutil.rmtree(tmp_csv_dir, ignore_errors=True)
            return
        env_pairs.append(f"{env_key}={val}")

    env_pairs += [
        f"METCELL_D_PATH=/data/{d_basename}",
        "METCELL_OUTPUT_DIR=/data",
        "METCELL_SCRIPT_DIR=/data",
        "METCELL_BRUKER_LIB=/data/libtimsdata.so",
        "METCELL_TIME_LIMIT_TABLE=/data/time_limit_table.csv",
    ]

    cmd_copy_in = [
        "docker", "run", "--rm",
        "-v", f"{d_parent_abs}:/src_d:ro",
        "-v", f"{metcell_abs}:/src_s:ro",
        "-v", f"{tmp_csv_dir}:/src_t:ro",
        "-v", f"{vol_name}:/dest",
        "ubuntu:24.04", "sh", "-c",
        f"cp -a /src_d/{d_basename} /dest/ && cp -a /src_s/. /dest/ "
        f"&& cp /src_t/time_limit_table.csv /dest/time_limit_table.csv"
    ]
    cmd_run = [
        "docker", "run", "--rm",
        "--name", container_name,
        "-w", "/data",
        "-v", f"{vol_name}:/data",
    ]
    for pair in env_pairs:
        cmd_run += ["-e", pair]
    cmd_run += [image, "Rscript", "/data/run.R"]

    cmd_copy_out = [
        "docker", "run", "--rm",
        "-v", f"{vol_name}:/data:ro",
        "-v", f"{out_abs}:/export",
        "ubuntu:24.04", "sh", "-c",
        "if [ -d /data/results ]; then cp -a /data/results/. /export/; fi"
    ]

    mc_result_text.delete(1.0, tk.END)
    mc_result_text.insert(tk.END,
        f"Using image : {image}\n"
        f"Volume      : {vol_name}\n"
        f"Container   : {container_name}\n"
        f"D folder    : {d_folder_abs}\n"
        f"Output dir  : {out_abs}\n"
        f"Segments    : {len(segments)}\n\n")
    root.update()

    btn_mc_run.config(state="disabled")
    if btn_worker_tims_run is not None:
        btn_worker_tims_run.config(state="disabled")
    btn_mc_terminate.config(state="normal")

    output_queue = queue.Queue()
    proc_holder = {'proc': None}

    def target():
        global metcell_proc

        def emit(msg):
            output_queue.put(("STDOUT", msg))

        try:
            emit(f"\n[0/4] Creating volume {vol_name}...")
            r = subprocess.run(["docker", "volume", "create", vol_name], capture_output=True,
                               text=True, encoding='utf-8', errors='replace', timeout=60)
            if r.returncode != 0:
                output_queue.put(("ERROR", f"volume create failed: {r.stderr}"))
                output_queue.put(("DONE", 1))
                return

            emit(f"\n[1/4] Copying data into volume...")
            emit("  cmd: " + " ".join(cmd_copy_in))
            rc = _docker_run_streaming(cmd_copy_in, output_queue, proc_holder)
            metcell_proc = proc_holder['proc']
            if rc != 0:
                output_queue.put(("ERROR", f"copy-in failed (rc={rc})"))
                output_queue.put(("DONE", rc))
                return

            emit(f"\n[2/4] Running MetCell...")
            emit("  cmd: " + " ".join(cmd_run))
            rc = _docker_run_streaming(cmd_run, output_queue, proc_holder)
            metcell_proc = proc_holder['proc']
            if rc != 0:
                output_queue.put(("ERROR", f"MetCell run failed (rc={rc})"))
                output_queue.put(("DONE", rc))
                return

            emit(f"\n[3/4] Copying results out...")
            emit("  cmd: " + " ".join(cmd_copy_out))
            rc = _docker_run_streaming(cmd_copy_out, output_queue, proc_holder)
            metcell_proc = proc_holder['proc']
            if rc != 0:
                output_queue.put(("ERROR", f"copy-out failed (rc={rc})"))
                output_queue.put(("DONE", rc))
                return

            emit(f"\n[4/4] Cleaning up volume {vol_name}...")
            _cleanup_volume(vol_name)
            output_queue.put(("DONE", 0))
        except Exception as e:
            output_queue.put(("ERROR", f"Exception: {e}"))
            output_queue.put(("DONE", 1))
        finally:
            shutil.rmtree(tmp_csv_dir, ignore_errors=True)

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    rc = item[1]
                    if rc == 0:
                        mc_result_text.insert(tk.END,
                            "\n*** MetCell completed successfully ***\n")
                        mc_result_text.insert(tk.END,
                            f"\nResults are in:\n  {os.path.join(out_abs, 'results')}\n")
                        btn_mc_run.config(state="normal")
                        btn_mc_terminate.config(state="disabled")
                        if btn_worker_tims_run is not None:
                            btn_worker_tims_run.config(state="normal")
                        # ✅ MetCell 成功后自动调用 worker_tims（其内部已集成 quantification）
                        root.after(100, run_worker_tims_subprocess)
                        return
                    else:
                        mc_result_text.insert(tk.END,
                            f"\n*** MetCell FAILED (rc={rc}) ***\n"
                            "Volume was NOT cleaned. You can inspect it with:\n"
                            f"  docker run --rm -v {vol_name}:/data ubuntu:24.04 ls /data\n")
                        btn_mc_run.config(state="normal")
                        btn_mc_terminate.config(state="disabled")
                        return
                elif item[0] == "ERROR":
                    mc_result_text.insert(tk.END, f"\n!!! {item[1]} !!!\n")
                    if btn_worker_tims_run is not None:
                            btn_worker_tims_run.config(state="normal")
                else:
                    prefix = ">> " if item[0] == "STDOUT" else "!! "
                    mc_result_text.insert(tk.END, f"{prefix}{item[1]}\n")
                    mc_result_text.see(tk.END)
        except queue.Empty:
            pass
        root.after(100, poll_queue)

    threading.Thread(target=target, daemon=True).start()
    root.after(100, poll_queue)


def terminate_metcell_subprocess():
    global metcell_proc, metcell_container_name, metcell_volume_name, worker_tims_proc
    any_running = (metcell_proc is not None) or (worker_tims_proc is not None)
    if not any_running:
        messagebox.showinfo("Info", "No MetCell / worker_tims process running")
        return
    try:
        if worker_tims_proc is not None:
            try:
                worker_tims_proc.terminate()
                try:
                    worker_tims_proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    worker_tims_proc.kill()
            except Exception:
                pass
            mc_result_text.insert(tk.END, "\n⚠️ User terminated worker_tims\n")

        if metcell_proc is not None:
            metcell_proc.terminate()
            try:
                metcell_proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                metcell_proc.kill()
            _cleanup_volume(metcell_volume_name, metcell_container_name)
            mc_result_text.insert(tk.END, "\n⚠️ User terminated MetCell (volume cleaned)\n")
    except Exception as e:
        messagebox.showerror("Termination Error", str(e))
    finally:
        metcell_proc = None
        metcell_container_name = None
        metcell_volume_name = None
        worker_tims_proc = None
        try:
            btn_mc_run.config(state="normal")
            if btn_worker_tims_run is not None:
                btn_worker_tims_run.config(state="normal")
            if btn_cs_tims_terminate is not None:
                btn_cs_tims_terminate.config(state="disabled")   # 新增
            btn_mc_terminate.config(state="disabled")
        except Exception:
            pass
def terminate_worker_tims_subprocess():
    """只终止 worker_tims（MPCell）进程，不影响 MetCell"""
    global worker_tims_proc
    if worker_tims_proc is None:
        messagebox.showinfo("Info", "No MPCell (worker_tims) process running")
        return
    try:
        worker_tims_proc.terminate()
        try:
            worker_tims_proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            worker_tims_proc.kill()
        mc_result_text.insert(tk.END, "\n⚠️ User terminated MPCell (worker_tims) only\n")
        mc_result_text.see(tk.END)
    except Exception as e:
        messagebox.showerror("Termination Error", str(e))
    finally:
        worker_tims_proc = None
        try:
            if btn_worker_tims_run is not None:
                btn_worker_tims_run.config(state="normal")
            if btn_cs_tims_terminate is not None:
                btn_cs_tims_terminate.config(state="disabled")
            # 若 MetCell 也不在跑，恢复 MetCell 页按钮；否则保持 disabled
            if metcell_proc is None:
                btn_mc_run.config(state="normal")
                btn_mc_terminate.config(state="disabled")
        except Exception:
            pass
def build_worker_tims_params():
    """从 mc_vars 和 mc_bg_rows 收集 worker_tims 需要的参数。
    d_folder 指向 MetCell 输出目录（因为 worker_tims 读的是 out_abs 下的 npy/csv）。
    """
    out_dir = mc_vars["output_dir"].get().strip()
    if not out_dir:
        return None
    out_abs = os.path.abspath(out_dir)

    bg_marker = []
    bg_intensity_max = []
    for row in mc_bg_rows:
        mz_str = row['mz'].get().strip()
        ccs_str = row['ccs'].get().strip()
        int_str = row['intensity'].get().strip()
        if not mz_str or not ccs_str:
            continue
        try:
            bg_marker.append([float(mz_str), float(ccs_str)])
            bg_intensity_max.append(float(int_str) if int_str else 0.0)
        except ValueError:
            continue

    def _f(key, default):
        try:
            return float(mc_vars[key].get())
        except Exception:
            return default

    def _i(key, default):
        try:
            return int(mc_vars[key].get())
        except Exception:
            return default
    # ---- 新增：评分参数（SNR / Match / RSD） ----
    def _b(key, default):
        try:
            return bool(mc_vars[key].get())
        except Exception:
            return default

    use_snr   = _b("use_snr",   True)
    use_match = _b("use_match", True)
    use_rsd   = _b("use_rsd",   True)

    snr_thres   = _f("snr_thres",   3.5)
    snr_lower   = _f("snr_lower",   2.0)
    match_thres = _f("match_thres", 0.67)
    match_lower = _f("match_lower", 0.3)
    rsd_thres   = _f("rsd_thres",   0.2)
    rsd_lower   = _f("rsd_lower",   0.1)
    ppm = _f("mz_tol_ppm", 20.0)
    mobility_range = _f("mobility_range", 0.05)
    marker_mz = _f("marker_mz", 760.5845)
    marker_ccs = _f("marker_mobility", 1.402)
    int_lower = _f("int_lower", 60000.0)
    cell_eic_span = _i("cell_eic_span", 5)
    cell_eic_span_lim = _i("cell_eic_span_lim", 1)
    spec_num_lim = _i("spec_num_lim", 1)
    rt_min = _f("rt_range_min", 0.0)
    rt_max = _f("rt_range_max", 1e9)
    cell_occur_thres = _f("cell_occurrence", 5.0) / 100.0
    filter_threshold = _f("threshold", 10.0)
    plot_option = mc_vars["plot_option"].get() if "plot_option" in mc_vars else "none"
    plot_rt_min = _f("plot_rt_min", 5.0)
    plot_rt_max = _f("plot_rt_max", 10.0)

    return {
        'd_folder': out_abs,                 # MetCell 输出目录（npy/csv 所在）
        'output_dir': out_abs,               # worker_tims 输出目录
        'mz_tol_ppm': ppm,
        'mobility_range': mobility_range,
        'ccs_range': mobility_range,
        'rt_range_min': rt_min,
        'rt_range_max': rt_max,
        'cell_eic_span': cell_eic_span,
        'cell_eic_span_lim': cell_eic_span_lim,
        'spec_num_lim': spec_num_lim,
        'cell_occurrence_threshold': cell_occur_thres,
        'marker_mz': marker_mz,
        'marker_ccs': marker_ccs,
        'int_lower': int_lower,
        'bg_marker': bg_marker,
        'bg_intensity_max': bg_intensity_max,
        'plot_option': plot_option,
        'plot_rt_min': plot_rt_min,
        'plot_rt_max': plot_rt_max,
        'filter_threshold': filter_threshold,
        'use_snr':     use_snr,
        'snr_thres':   snr_thres,
        'snr_lower':   snr_lower,
        'use_match':   use_match,
        'match_thres': match_thres,
        'match_lower': match_lower,
        'use_rsd':     use_rsd,
        'rsd_thres':   rsd_thres,
        'rsd_lower':   rsd_lower,
    }


def run_worker_tims_subprocess():
    """MetCell 成功后调用 worker_tims。复用 mc_result_text 显示输出。"""
    global worker_tims_proc
    params = build_worker_tims_params()
    if params is None:
        mc_result_text.insert(tk.END, "\n[worker_tims] ERROR: output_dir not set\n")
        return

    out_abs = params['output_dir']
    param_file_path = os.path.join(out_abs, 'temp_worker_tims_params.json')
    try:
        with open(param_file_path, 'w', encoding='utf-8') as f:
            json.dump(params, f)
    except Exception as e:
        mc_result_text.insert(tk.END, f"\n[worker_tims] Cannot write params: {e}\n")
        return

    mc_result_text.insert(tk.END, "\n========== Starting worker_tims ==========\n")
    mc_result_text.see(tk.END)
    root.update()

    try:
        if btn_worker_tims_run is not None:
            btn_worker_tims_run.config(state="disabled")
        if btn_cs_tims_terminate is not None:
            btn_cs_tims_terminate.config(state="normal")         # 新增
        btn_mc_run.config(state="disabled")
        btn_mc_terminate.config(state="normal")
    except Exception:
        pass

    output_queue = queue.Queue()

    def read_stream(stream, name):
        try:
            for line in iter(stream.readline, ''):
                if line:
                    output_queue.put((name, line.rstrip()))
                else:
                    break
        except Exception as e:
            output_queue.put((name, f"Read error: {e}"))
        finally:
            stream.close()

    def target():
        global worker_tims_proc
        try:
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            proc = subprocess.Popen(
                get_worker_tims_command() + [param_file_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, bufsize=1, env=env
            )
            worker_tims_proc = proc
            t_out = threading.Thread(target=read_stream, args=(proc.stdout, "STDOUT"))
            t_err = threading.Thread(target=read_stream, args=(proc.stderr, "STDERR"))
            t_out.daemon = t_err.daemon = True
            t_out.start(); t_err.start()
            rc = proc.wait()
            t_out.join(timeout=1); t_err.join(timeout=1)
            output_queue.put(("DONE", rc))
        except Exception as e:
            output_queue.put(("ERROR", str(e)))

    def poll_queue():
        try:
            while True:
                item = output_queue.get_nowait()
                if item[0] == "DONE":
                    rc = item[1]
                    if rc == 0:
                        mc_result_text.insert(tk.END,
                            "\n*** worker_tims completed successfully ***\n")
                    else:
                        mc_result_text.insert(tk.END,
                            f"\n*** worker_tims FAILED (rc={rc}) ***\n")
                    try:
                        os.remove(param_file_path)
                    except Exception:
                        pass
                    try:
                        if btn_worker_tims_run is not None:
                            btn_worker_tims_run.config(state="normal")
                        if btn_cs_tims_terminate is not None:
                            btn_cs_tims_terminate.config(state="disabled")
                        btn_mc_run.config(state="normal")
                        btn_mc_terminate.config(state="disabled")
                    except Exception:
                        pass
                    mc_result_text.see(tk.END)
                    return
                elif item[0] == "ERROR":
                    mc_result_text.insert(tk.END, f"\n!!! worker_tims error: {item[1]} !!!\n")
                else:
                    prefix = ">> " if item[0] == "STDOUT" else "!! "
                    mc_result_text.insert(tk.END, f"[worker_tims] {prefix}{item[1]}\n")
                    mc_result_text.see(tk.END)
        except queue.Empty:
            pass
        root.after(100, poll_queue)

    threading.Thread(target=target, daemon=True).start()
    root.after(100, poll_queue)

# ============================================================
# ===============  构建 UI  =================================
# ============================================================
root = tk.Tk()
root.title("MPCell")
root.geometry("1200x850")
root.minsize(500, 500)
root.resizable(True, True)

status_label = ttk.Label(root, text="Ready", relief="sunken", anchor="w")
status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=2)

notebook = ttk.Notebook(root)
notebook.pack(fill='both', expand=True, padx=10, pady=10)

# ---------- 页 1：File Settings（始终显示）----------
# ---------- 页 1：File Settings（始终显示，仅保留 MS 模式选择器）----------
file_settings_sf = ScrollableFrame(notebook)
frame_files = file_settings_sf.inner

row = 0
ttk.Label(frame_files, text="MS:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
mode_var = tk.StringVar(value="no_tims")
frame_mode_sel = ttk.Frame(frame_files)
frame_mode_sel.grid(row=row, column=1, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_mode_sel, text="no_IM", variable=mode_var, value="no_tims").pack(side='left', padx=2)
ttk.Radiobutton(frame_mode_sel, text="IM-based", variable=mode_var, value="tims").pack(side='left', padx=2)
row += 1

frame_files.columnconfigure(1, weight=1)


# ---------- 页 2：Parameter Setting（no_tims 模式）----------
param_setting_sf = ScrollableFrame(notebook)
frame_basic = param_setting_sf.inner

row = 0

# ---- 文件路径（从 File Settings 页移来）----
ttk.Label(frame_basic, text="Input file (.mzML):").grid(
    row=row, column=0, sticky='e', padx=(0, 2), pady=5)
entry_input_file = ttk.Entry(frame_basic, width=80)
entry_input_file.grid(row=row, column=1, padx=(0, 5), pady=5, sticky='w')
ttk.Button(frame_basic, text="Browse",
           command=lambda: browse_input_file(entry_input_file)).grid(
    row=row, column=2, padx=5, pady=5)
row += 1

ttk.Label(frame_basic, text="Output direction:").grid(
    row=row, column=0, sticky='e', padx=(0, 2), pady=5)
entry_output_dir = ttk.Entry(frame_basic, width=80)
entry_output_dir.grid(row=row, column=1, padx=(0, 5), pady=5, sticky='w')
ttk.Button(frame_basic, text="Browse",
           command=lambda: browse_output_dir(entry_output_dir)).grid(
    row=row, column=2, padx=5, pady=5)
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

ttk.Label(frame_basic, text="【Basic setting】").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
row += 1

ttk.Label(frame_basic, text="Ion mode:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
frame_mode_ion = ttk.Frame(frame_basic)
frame_mode_ion.grid(row=row, column=1, columnspan=2, sticky='w', padx=0, pady=5)
ion_mode_var = tk.StringVar(value="Positive")
ttk.Radiobutton(frame_mode_ion, text="Positive", variable=ion_mode_var, value="Positive").pack(side='left', padx=5)
ttk.Radiobutton(frame_mode_ion, text="Negative", variable=ion_mode_var, value="Negative").pack(side='left', padx=5)
row += 1

ttk.Label(frame_basic, text="Peak intensity threshold:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_peak_intensity = ttk.Entry(frame_basic, width=20)
entry_peak_intensity.grid(row=row, column=1, padx=(0, 5), pady=3, sticky='w')
entry_peak_intensity.insert(0, "1e4")
row += 1

ttk.Label(frame_basic, text="mz tolerance [ppm]:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_ppm = ttk.Entry(frame_basic, width=20)
entry_ppm.grid(row=row, column=1, padx=(0, 5), pady=3, sticky='w')
entry_ppm.insert(0, "10")
row += 1

ttk.Label(frame_basic, text="RT range:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
frame_rt = ttk.Frame(frame_basic)
frame_rt.grid(row=row, column=1, sticky='w', padx=0, pady=3)
entry_rt_min = ttk.Entry(frame_rt, width=8); entry_rt_min.pack(side='left'); entry_rt_min.insert(0, "0")
ttk.Label(frame_rt, text="~").pack(side='left', padx=2)
entry_rt_max = ttk.Entry(frame_rt, width=8); entry_rt_max.pack(side='left'); entry_rt_max.insert(0, "100")
ttk.Label(frame_rt, text="min").pack(side='left', padx=2)
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

ttk.Label(frame_basic, text="【Cell marker setting】").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
row += 1

ttk.Label(frame_basic, text="Marker m/z:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_marker_mz = ttk.Entry(frame_basic, width=15)
entry_marker_mz.grid(row=row, column=1, sticky='w', padx=(0, 5), pady=3)
row += 1

ttk.Label(frame_basic, text="Marker intensity threshold:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_intensity_thres = ttk.Entry(frame_basic, width=20)
entry_intensity_thres.grid(row=row, column=1, sticky='w', padx=(0, 5), pady=3)
row += 1

ttk.Label(frame_basic, text="Cell EIC peak span:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_cell_eic_span = ttk.Entry(frame_basic, width=20)
entry_cell_eic_span.grid(row=row, column=1, padx=(0, 5), pady=3, sticky='w')
entry_cell_eic_span.insert(0, "5")
ttk.Label(frame_basic, text="(default 5)", foreground="gray").grid(row=row, column=2, sticky='w', padx=2, pady=3)
row += 1

ttk.Label(frame_basic, text="Limited Cell EIC peak span:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
entry_cell_eic_span_lim = ttk.Entry(frame_basic, width=20)
entry_cell_eic_span_lim.grid(row=row, column=1, padx=(0, 5), pady=3, sticky='w')
entry_cell_eic_span_lim.insert(0, "1")
ttk.Label(frame_basic, text="(default 1, recommand >=2 when cell pulse spectra number > 10)",
          foreground="gray").grid(row=row, column=2, sticky='w', padx=2, pady=3)
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

ttk.Label(frame_basic, text="【Background marker setting】").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
row += 1
ttk.Label(frame_basic, text="Input format: m/z,max intensity").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
row += 1

bg_rows_container = tk.Frame(frame_basic)
bg_rows_container.grid(row=row, column=1, sticky='ew', pady=5)

btn_frame_bg = ttk.Frame(frame_basic)
btn_frame_bg.grid(row=row, column=2, padx=5, sticky='n')
ttk.Button(btn_frame_bg, text="Add row", command=lambda: add_bg_row()).pack(pady=2)
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

ttk.Label(frame_basic, text="【EIC plot setting】").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
row += 1

ttk.Label(frame_basic, text="EIC Plot:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
plot_var = tk.StringVar(value="none")
frame_plot = ttk.Frame(frame_basic)
frame_plot.grid(row=row, column=1, columnspan=2, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_plot, text="None", variable=plot_var, value="none").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="Cell-specific peak only", variable=plot_var, value="fp").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="Non-specific peak only", variable=plot_var, value="trash").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot, text="All", variable=plot_var, value="all").pack(side='left', padx=2)
row += 1

ttk.Label(frame_basic, text="Plot RT range:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=3)
frame_rt_plot = ttk.Frame(frame_basic)
frame_rt_plot.grid(row=row, column=1, sticky='w', padx=0, pady=3)
entry_plot_rt_min = ttk.Entry(frame_rt_plot, width=8); entry_plot_rt_min.pack(side='left'); entry_plot_rt_min.insert(0, "5")
ttk.Label(frame_rt_plot, text="~").pack(side='left', padx=2)
entry_plot_rt_max = ttk.Entry(frame_rt_plot, width=8); entry_plot_rt_max.pack(side='left'); entry_plot_rt_max.insert(0, "10")
ttk.Label(frame_rt_plot, text="min").pack(side='left', padx=2)
row += 1

ttk.Separator(frame_basic, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky='ew', pady=8)
row += 1

frame_basic.columnconfigure(1, weight=1)

# ---------- 页 3：Cell-specific Peak Parameters ----------
cell_specific_sf = ScrollableFrame(notebook)
frame_cs = cell_specific_sf.inner

row = 0
ttk.Label(frame_cs, text="【Cell-specific peak filtration】",
          font=('Arial', 11, 'bold')).grid(
    row=row, column=0, columnspan=5, sticky='w', padx=5, pady=(10, 6))
row += 1

# 表头
ttk.Label(frame_cs, text="Use",       font=('Arial', 10, 'bold')).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs, text="Parameter", font=('Arial', 10, 'bold')).grid(row=row, column=1, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs, text="Threshold", font=('Arial', 10, 'bold')).grid(row=row, column=2, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs, text="Lower",     font=('Arial', 10, 'bold')).grid(row=row, column=3, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs, text="Default",   font=('Arial', 10, 'bold')).grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# SNR
var_use_snr = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs, variable=var_use_snr).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs, text="SNR:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
entry_snr_thres = ttk.Entry(frame_cs, width=10)
entry_snr_thres.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_snr_thres.insert(0, "3.5")
entry_snr_lower = ttk.Entry(frame_cs, width=10)
entry_snr_lower.grid(row=row, column=3, sticky='w', padx=8, pady=3)
entry_snr_lower.insert(0, "2")
ttk.Label(frame_cs, text="(3.5 / 2)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# Match ratio
var_use_match = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs, variable=var_use_match).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs, text="Match ratio:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
entry_match_thres = ttk.Entry(frame_cs, width=10)
entry_match_thres.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_match_thres.insert(0, "0.67")
entry_match_lower = ttk.Entry(frame_cs, width=10)
entry_match_lower.grid(row=row, column=3, sticky='w', padx=8, pady=3)
entry_match_lower.insert(0, "0.3")
ttk.Label(frame_cs, text="(0.67 / 0.3)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# RSD subs
var_use_rsd = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs, variable=var_use_rsd).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs, text="RSD subs:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
entry_rsd_thres = ttk.Entry(frame_cs, width=10)
entry_rsd_thres.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_rsd_thres.insert(0, "0.2")
entry_rsd_lower = ttk.Entry(frame_cs, width=10)
entry_rsd_lower.grid(row=row, column=3, sticky='w', padx=8, pady=3)
entry_rsd_lower.insert(0, "0.1")
ttk.Label(frame_cs, text="(0.2 / 0.1)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

ttk.Separator(frame_cs, orient='horizontal').grid(
    row=row, column=0, columnspan=5, sticky='ew', pady=10)
row += 1

# Score threshold
ttk.Label(frame_cs, text="Score threshold:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
entry_threshold = ttk.Entry(frame_cs, width=10)
entry_threshold.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_threshold.insert(0, "10")
threshold_hint = ttk.Label(frame_cs, text="(default 10)", foreground="gray")
threshold_hint.grid(row=row, column=3, sticky='w', padx=8, pady=3)
entry_threshold.bind('<KeyRelease>', on_threshold_edit)
threshold_is_default = True
row += 1

# Limitation of detected spectra in cell pulse
ttk.Label(frame_cs, text="Limitation of detected spectra in cell pulse:").grid(
    row=row, column=1, sticky='w', padx=8, pady=3)
entry_spec_num_lim = ttk.Entry(frame_cs, width=10)
entry_spec_num_lim.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_spec_num_lim.insert(0, "1")
ttk.Label(frame_cs, text="(default 1, recommand 2 when cell pulse spectra number > 10)",
          foreground="gray").grid(row=row, column=3, sticky='w', padx=8, pady=3)
row += 1

# Cell occurrence
ttk.Label(frame_cs, text="Cell occurrence (%):").grid(row=row, column=1, sticky='w', padx=8, pady=3)
entry_cell_occurrence = ttk.Entry(frame_cs, width=10)
entry_cell_occurrence.grid(row=row, column=2, sticky='w', padx=8, pady=3)
entry_cell_occurrence.insert(0, "5")
ttk.Label(frame_cs, text="(default 5)", foreground="gray").grid(row=row, column=3, sticky='w', padx=8, pady=3)
row += 1

ttk.Separator(frame_cs, orient='horizontal').grid(
    row=row, column=0, columnspan=5, sticky='ew', pady=10)
row += 1

# Run / Terminate 按钮（放本页）
btn_frame_run = ttk.Frame(frame_cs)
btn_frame_run.grid(row=row, column=0, columnspan=5, pady=15)
btn_run = ttk.Button(btn_frame_run, text="Run",
                     command=run_subprocess_threaded, width=15)
btn_run.pack(side='left', padx=10)
btn_terminate = ttk.Button(btn_frame_run, text="Terminate",
                           command=terminate_subprocess, width=15, state="disabled")
btn_terminate.pack(side='left', padx=10)
row += 1

# ---------- 页 3：Output（no_tims 模式）----------
output_sf = ScrollableFrame(notebook)
frame_output = output_sf.inner

result_text = tk.Text(frame_output, wrap=tk.WORD, height=60)
scrollbar_result = ttk.Scrollbar(frame_output, command=result_text.yview)
result_text.configure(yscrollcommand=scrollbar_result.set)
result_text.pack(side='left', fill='both', expand=True, padx=5, pady=5)
scrollbar_result.pack(side='right', fill='y')


# ---------- 页 4：Quantification（no_tims 模式）----------
quant_sf = ScrollableFrame(notebook)
frame_quant = quant_sf.inner

row = 0
ttk.Label(frame_quant, text="Select peak_list file (.xlsx):").grid(row=row, column=0, sticky='w', padx=5, pady=5)
entry_mz_list = ttk.Entry(frame_quant, width=60)
entry_mz_list.grid(row=row, column=1, padx=(0, 5), pady=5, sticky='w')
ttk.Button(frame_quant, text="Browse",
           command=lambda: browse_mz_list(entry_mz_list)).grid(row=row, column=2, padx=5, pady=5)
row += 1

ttk.Label(frame_quant, text="Note: Excel should contain a 'mz' (or 'm/z') column").grid(
    row=row, column=0, sticky='w', padx=5, pady=2)
row += 1
ttk.Label(frame_quant, text="(default dir: /grouping/RT_del/cell_specific_peak.xlsx)",
          foreground="gray").grid(row=row, column=0, sticky='w', padx=5, pady=5)
row += 1

btn_frame_quant = ttk.Frame(frame_quant)
btn_frame_quant.grid(row=row, column=0, columnspan=3, pady=10)
btn_quant_run = ttk.Button(btn_frame_quant, text="Run", command=run_quant_subprocess, width=12)
btn_quant_run.pack(side='left', padx=5)
btn_quant_terminate = ttk.Button(btn_frame_quant, text="Terminate",
                                 command=terminate_quant_subprocess, width=12, state="disabled")
btn_quant_terminate.pack(side='left', padx=5)
row += 1

ttk.Label(frame_quant, text="Quantification Output:").grid(row=row, column=0, sticky='w', padx=5, pady=2)
row += 1
quant_result_text = tk.Text(frame_quant, wrap=tk.WORD, height=20)
scrollbar_quant = ttk.Scrollbar(frame_quant, command=quant_result_text.yview)
quant_result_text.configure(yscrollcommand=scrollbar_quant.set)
quant_result_text.grid(row=row, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)
scrollbar_quant.grid(row=row, column=3, sticky='ns')
frame_quant.columnconfigure(1, weight=1)


# ---------- 页 5：Metabolite Annotation（no_tims 模式）----------
anno_sf = ScrollableFrame(notebook)
frame_annotation = anno_sf.inner
frame_annotation.columnconfigure(0, weight=1)
frame_annotation.columnconfigure(1, weight=1)

row = 0
ttk.Label(frame_annotation, text="Ion mode:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
frame_basic_anno = ttk.Frame(frame_annotation)
frame_basic_anno.grid(row=row, column=1, sticky='w', pady=5, columnspan=1)
ion_mode_var_anno = tk.StringVar(value="Positive")
ttk.Radiobutton(frame_basic_anno, text="Positive", variable=ion_mode_var_anno, value="Positive").pack(side='left', padx=2)
ttk.Radiobutton(frame_basic_anno, text="Negative", variable=ion_mode_var_anno, value="Negative").pack(side='left', padx=2)
row += 1

ttk.Label(frame_annotation, text="Input m/z file (.xlsx):").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
entry_mz_file = ttk.Entry(frame_annotation, width=60)
entry_mz_file.grid(row=row, column=1, padx=(0, 5), pady=5, sticky='w')
ttk.Button(frame_annotation, text="Browse",
           command=lambda: browse_mz_file(entry_mz_file)).grid(row=row, column=2, padx=5, pady=5)
row += 1

ttk.Label(frame_annotation, text="(default dir: /intensity_matrix_iso_anno_M0.xlsx)",
          foreground="gray").grid(row=row, column=1, sticky='w', padx=5, pady=5)
row += 1

ttk.Label(frame_annotation, text="Embedded databases:").grid(row=row, column=0, sticky='ne', padx=(0, 2), pady=5)
frame_embedded = ttk.Frame(frame_annotation)
frame_embedded.grid(row=row, column=1, columnspan=2, sticky='w', padx=(0, 2), pady=5)
embedded_names = ["lcms_mda", "lcms_hela", "lcms_ninecells", "sc_paper", "hmdb", "lipidmaps"]
for i, name in enumerate(embedded_names):
    var = tk.BooleanVar(value=True)
    embedded_db_vars[name] = var
    ttk.Checkbutton(frame_embedded, text=name, variable=var,
                    command=update_db_priority_list).grid(row=i // 3, column=i % 3, sticky='w', padx=5, pady=2)
row += 1

ttk.Label(frame_annotation, text="Custom database:").grid(row=row, column=0, sticky='e', padx=(0, 2), pady=5)
entry_db_file = ttk.Entry(frame_annotation, width=50)
entry_db_file.grid(row=row, column=1, padx=(0, 5), pady=5, sticky='w')
ttk.Button(frame_annotation, text="Add", command=add_custom_db, width=6).grid(row=row, column=2, padx=2, pady=5, sticky='w')
ttk.Button(frame_annotation, text="Remove", command=remove_custom_db, width=10).grid(row=row, column=3, padx=2, pady=5, sticky='w')
row += 1

priority_frame = ttk.Frame(frame_annotation)
priority_frame.grid(row=row, column=0, columnspan=4, sticky='nsew', pady=5)
priority_frame.columnconfigure(0, weight=1)
priority_frame.columnconfigure(1, weight=1)

left_frame = ttk.Frame(priority_frame, relief='ridge', borderwidth=1)
left_frame.grid(row=0, column=0, sticky='nsew', padx=5)
ttk.Label(left_frame, text="Database priority", font=('Arial', 11, 'bold')).pack(pady=2)

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

db_priority_list = [name for name in embedded_names]
refresh_db_listbox()

right_frame = ttk.Frame(priority_frame, relief='ridge', borderwidth=1)
right_frame.grid(row=0, column=1, sticky='nsew', padx=5)
ttk.Label(right_frame, text="Adduct priority", font=('Arial', 11, 'bold')).pack(pady=2)

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

ion_mode_var_anno.trace('w', on_ion_mode_change)
init_adduct_priority(ion_mode_var_anno.get())
row += 1

param_frame_anno = ttk.Frame(frame_annotation)
param_frame_anno.grid(row=row, column=0, columnspan=4, sticky='w', padx=5, pady=5)

ttk.Label(param_frame_anno, text="m/z tolerance:").pack(side='left', padx=(0, 2))
entry_ppm_anno = ttk.Entry(param_frame_anno, width=10)
entry_ppm_anno.pack(side='left', padx=(0, 10))
entry_ppm_anno.insert(0, "10")

ttk.Label(param_frame_anno, text="Max candidates:").pack(side='left', padx=(0, 2))
entry_max_candidates = ttk.Entry(param_frame_anno, width=10)
entry_max_candidates.pack(side='left')
entry_max_candidates.insert(0, "5")
row += 1

btn_frame_anno = ttk.Frame(frame_annotation)
btn_frame_anno.grid(row=row, column=0, columnspan=4, pady=10)
btn_anno_run = ttk.Button(btn_frame_anno, text="Run", command=run_annotation_subprocess, width=12)
btn_anno_run.pack(side='left', padx=5)
btn_anno_terminate = ttk.Button(btn_frame_anno, text="Terminate",
                                command=terminate_annotation_subprocess, width=12, state="disabled")
btn_anno_terminate.pack(side='left', padx=5)
row += 1

ttk.Label(frame_annotation, text="Annotation Output:").grid(row=row, column=0, sticky='w', padx=5, pady=2)
row += 1
anno_result_text = tk.Text(frame_annotation, wrap=tk.WORD, height=10)
scrollbar_anno = ttk.Scrollbar(frame_annotation, command=anno_result_text.yview)
anno_result_text.configure(yscrollcommand=scrollbar_anno.set)
anno_result_text.grid(row=row, column=0, columnspan=4, sticky='nsew', padx=5, pady=5)
scrollbar_anno.grid(row=row, column=4, sticky='ns')


# ---------- 页 6：MetCell EIC（tims 模式，所有变量在 mc_vars 中）----------
metcell_sf = ScrollableFrame(notebook)
frame_metcell = metcell_sf.inner

frame_metcell.columnconfigure(0, weight=0, minsize=130)
frame_metcell.columnconfigure(1, weight=1)
frame_metcell.columnconfigure(2, weight=0)

mc_row = 0

ttk.Label(frame_metcell, text="【Paths】", font=('Arial', 11, 'bold')).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(10, 4))
mc_row += 1


def _mc_add_path(label, var_key, default="", mode="dir"):
    global mc_row
    ttk.Label(frame_metcell, text=label).grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
    mc_vars[var_key] = tk.StringVar(value=default)
    ttk.Entry(frame_metcell, width=60, textvariable=mc_vars[var_key]).grid(
        row=mc_row, column=1, padx=(0, 5), pady=3, sticky='w')
    if mode == "dir":
        cmd = lambda: _browse_into(mc_vars[var_key], True)
    elif mode == "file":
        cmd = lambda: _browse_into(mc_vars[var_key], False)
    else:
        cmd = None
    if cmd:
        ttk.Button(frame_metcell, text="Browse", command=cmd).grid(row=mc_row, column=2, padx=5, pady=3)
    mc_row += 1


_mc_add_path("Bruker .d folder:", "d_folder", mode="dir")
_mc_add_path("Output directory:", "output_dir", mode="dir")

ttk.Label(frame_metcell, text="【Time segments (time_limit_table.csv)】",
          font=('Arial', 11, 'bold')).grid(row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(12, 4))
mc_row += 1

mc_vars["data_file"] = tk.StringVar(value="")
ttk.Label(frame_metcell, text="data_file (auto):").grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
ttk.Entry(frame_metcell, width=40, textvariable=mc_vars["data_file"], state='readonly').grid(
    row=mc_row, column=1, padx=(0, 5), pady=3, sticky='w')
mc_row += 1

seg_outer_frame = ttk.Frame(frame_metcell)
seg_outer_frame.grid(row=mc_row, column=0, columnspan=3, sticky='ew', padx=5, pady=3)
mc_row += 1

seg_header = ttk.Frame(seg_outer_frame)
seg_header.pack(fill='x')
ttk.Label(seg_header, text="start_time (s)", width=15).pack(side='left', padx=2)
ttk.Label(seg_header, text="end_time (s)", width=15).pack(side='left', padx=2)

seg_rows_container = tk.Frame(seg_outer_frame)
seg_rows_container.pack(fill='x', pady=2)

seg_btn_bar = ttk.Frame(seg_outer_frame)
seg_btn_bar.pack(fill='x', pady=2)
ttk.Button(seg_btn_bar, text="Add segment",
           command=lambda: add_time_seg_row()).pack(side='left', padx=2)
ttk.Button(seg_btn_bar, text="Clear all", command=clear_time_seg_rows).pack(side='left', padx=2)

mc_vars["d_folder"].trace('w', lambda *a: update_data_file_display())

ttk.Label(frame_metcell, text="Docker image:").grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
mc_vars["docker_image"] = tk.StringVar(value="zhulab/metcell-r:latest")
ttk.Entry(frame_metcell, width=60, textvariable=mc_vars["docker_image"]).grid(
    row=mc_row, column=1, padx=(0, 5), pady=3, sticky='w')
mc_row += 1

ttk.Label(frame_metcell, text="【Basic parameters】", font=('Arial', 11, 'bold')).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(12, 4))
mc_row += 1

ttk.Label(frame_metcell, text="Ion mode:").grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
mc_vars["ion_mode"] = tk.StringVar(value="positive")
_frame_basic_mc = ttk.Frame(frame_metcell)
_frame_basic_mc.grid(row=mc_row, column=1, sticky='w', pady=3)
ttk.Radiobutton(_frame_basic_mc, text="positive",
                variable=mc_vars["ion_mode"], value="positive").pack(side='left', padx=2)
ttk.Radiobutton(_frame_basic_mc, text="negative",
                variable=mc_vars["ion_mode"], value="negative").pack(side='left', padx=2)
mc_row += 1


def _mc_add_param(label, var_key, default):
    global mc_row
    ttk.Label(frame_metcell, text=label).grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
    mc_vars[var_key] = tk.StringVar(value=default)
    ttk.Entry(frame_metcell, width=20, textvariable=mc_vars[var_key]).grid(
        row=mc_row, column=1, padx=(0, 5), pady=3, sticky='w')
    mc_row += 1


_mc_add_param("nSlaves (threads):", "n_slaves", "2")
ttk.Label(frame_metcell, text="RT range (s):").grid(
    row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
frame_rt_range_mc = ttk.Frame(frame_metcell)
frame_rt_range_mc.grid(row=mc_row, column=1, sticky='w', pady=3)
mc_vars["rt_range_min"] = tk.StringVar(value="0")
mc_vars["rt_range_max"] = tk.StringVar(value="7200")
ttk.Entry(frame_rt_range_mc, width=8,
          textvariable=mc_vars["rt_range_min"]).pack(side='left')
ttk.Label(frame_rt_range_mc, text="~").pack(side='left', padx=2)
ttk.Entry(frame_rt_range_mc, width=8,
          textvariable=mc_vars["rt_range_max"]).pack(side='left')
ttk.Label(frame_rt_range_mc, text="s").pack(side='left', padx=2)
mc_row += 1

ttk.Label(frame_metcell, text="【Marker / Single-cell detection】", font=('Arial', 11, 'bold')).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(12, 4))
mc_row += 1

_mc_add_param("Marker m/z:", "marker_mz", "760.5845")
_mc_add_param("Marker mobility:", "marker_mobility", "1.402")
_mc_add_param("m/z tolerance [ppm]:", "mz_tol_ppm", "20")
_mc_add_param("Mobility range:", "mobility_range", "0.05")
_mc_add_param("Intensity upper:", "int_upper", "600000")
_mc_add_param("Intensity lower:", "int_lower", "60000")
_mc_add_param("Marker EIM peak span:", "marker_eim_span", "5")
_mc_add_param("Marker EIC peak span:", "cell_eic_span", "1")
_mc_add_param("Limited marker EIC peak span:", "cell_eic_span_lim", "1")

ttk.Separator(frame_metcell, orient='horizontal').grid(row=mc_row, column=0, columnspan=2, sticky='ew', pady=8)
mc_row += 1

ttk.Label(frame_metcell, text="【Background marker setting】", font=('Arial', 11, 'bold')).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(10, 4))
mc_row += 1

ttk.Label(frame_metcell, text="Input format: m/z,mobility,max intensity").grid(
    row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
mc_row += 1

mc_bg_rows_container = ttk.Frame(frame_metcell)
mc_bg_rows_container.grid(row=mc_row, column=1, sticky='ew', pady=5)

btn_frame_bg_mc = ttk.Frame(frame_metcell)
btn_frame_bg_mc.grid(row=mc_row, column=2, padx=5, sticky='n')
ttk.Button(btn_frame_bg_mc, text="Add row",
           command=lambda: add_mc_bg_row()).pack(pady=2)
mc_row += 1

ttk.Separator(frame_metcell, orient='horizontal').grid(row=mc_row, column=0, columnspan=2, sticky='ew', pady=8)
mc_row += 1

ttk.Label(frame_metcell, text="【EIC plot setting】", font=('Arial', 11, 'bold')).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(10, 4))
mc_row += 1

ttk.Label(frame_metcell, text="EIC Plot:").grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=5)
mc_vars["plot_option"] = tk.StringVar(value="none")
frame_plot_mc = ttk.Frame(frame_metcell)
frame_plot_mc.grid(row=mc_row, column=1, columnspan=2, sticky='w', padx=0, pady=5)
ttk.Radiobutton(frame_plot_mc, text="None", variable=mc_vars["plot_option"], value="none").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot_mc, text="Cell-specific peak only",
                variable=mc_vars["plot_option"], value="fp").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot_mc, text="Non-specific peak only",
                variable=mc_vars["plot_option"], value="trash").pack(side='left', padx=2)
ttk.Radiobutton(frame_plot_mc, text="All", variable=mc_vars["plot_option"], value="all").pack(side='left', padx=2)
mc_row += 1

ttk.Label(frame_metcell, text="Plot RT range:").grid(row=mc_row, column=0, sticky='e', padx=(0, 2), pady=3)
frame_rt_plot_mc = ttk.Frame(frame_metcell)
frame_rt_plot_mc.grid(row=mc_row, column=1, sticky='w', pady=3)
mc_vars["plot_rt_min"] = tk.StringVar(value="5")
mc_vars["plot_rt_max"] = tk.StringVar(value="10")
ttk.Entry(frame_rt_plot_mc, width=8, textvariable=mc_vars["plot_rt_min"]).pack(side='left')
ttk.Label(frame_rt_plot_mc, text="~").pack(side='left', padx=2)
ttk.Entry(frame_rt_plot_mc, width=8, textvariable=mc_vars["plot_rt_max"]).pack(side='left')
ttk.Label(frame_rt_plot_mc, text="min").pack(side='left', padx=2)
mc_row += 1

ttk.Separator(frame_metcell, orient='horizontal').grid(row=mc_row, column=0, columnspan=2, sticky='ew', pady=8)
mc_row += 1

# 高级参数（折叠）
mc_adv_visible = tk.BooleanVar(value=False)


def _toggle_mc_adv():
    if mc_adv_visible.get():
        frame_mc_adv.grid()
    else:
        frame_mc_adv.grid_remove()


ttk.Checkbutton(frame_metcell, text="Show advanced parameters",
                variable=mc_adv_visible, command=_toggle_mc_adv).grid(
    row=mc_row, column=0, columnspan=3, sticky='w', padx=5, pady=(10, 4))
mc_row += 1

frame_mc_adv = ttk.Frame(frame_metcell)
frame_mc_adv.grid(row=mc_row, column=0, columnspan=3, sticky='ew', padx=5, pady=2)
mc_row += 1
frame_mc_adv.grid_remove()

ADV_PARAMS = [
    ("Query pool size", "query_pool", "2000"),
    ("Superposition m/z tol [ppm]", "super_mz_tol", "20"),
    ("Cell number", "cell_number", "100"),
    ("Superposition pool", "super_pool", "40"),
    ("Search m/z tol [ppm]", "search_mz_tol", "20"),
    ("Search n_skip", "search_n_skip", "1"),
    ("Peak target length", "peak_len", "15"),
    ("Smooth window", "smooth_window", "10"),
    ("Peak span detection", "span_det", "21"),
    ("Peak span integration", "span_intg", "27"),
    ("Signal SD threshold", "signal_sd", "0.35"),
    ("Single-charge slope", "chg_slope", "0.0009"),
    ("Single-charge intercept", "chg_intercept", "0.35"),
    ("Extract m/z tol [ppm]", "extract_mz", "20"),
    ("Mobility tol", "extract_mob", "0.04"),
    ("Min cell number", "min_cell_num", "10"),
    ("Extract pool", "extract_pool", "150"),
    ("Isotope m/z tol [ppm]", "iso_mz_tol", "20"),
    ("Isotope mobility tol", "iso_mob_tol", "0.01"),
    ("Isotope delta [Da]", "iso_delta", "1.003355"),
    ("Isotope max num", "iso_max_num", "4"),
    ("Isotope ratio check (T/F)", "iso_ratio_chk", "TRUE"),
    ("Isotope ratio cutoff", "iso_ratio_cut", "500"),
    ("SC type CCS", "sc_type_ccs", "percentage"),
    ("SC CCS tol", "sc_tol_ccs", "4"),
    ("SC m/z tol [ppm]", "sc_tol_mz", "20"),
    ("Adduct for id", "sc_adduct", "[M+H]+"),
    ("EIC m/z tol [ppm]", "eic_mz_tol", "20"),
    ("EIC chunk size", "eic_chunk", "500"),
    ("EIC n threads", "eic_n_threads", "4"),
]

for i, (lbl, key, dflt) in enumerate(ADV_PARAMS):
    r = i // 3
    c = (i % 3) * 2
    ttk.Label(frame_mc_adv, text=lbl + ":").grid(row=r, column=c, sticky='e', padx=(4, 2), pady=2)
    mc_vars[key] = tk.StringVar(value=dflt)
    ttk.Entry(frame_mc_adv, width=10, textvariable=mc_vars[key]).grid(
        row=r, column=c + 1, sticky='w', padx=(0, 6), pady=2)

btn_frame_mc = ttk.Frame(frame_metcell)
btn_frame_mc.grid(row=mc_row, column=0, columnspan=3, pady=10)
mc_row += 1
btn_mc_run = ttk.Button(btn_frame_mc, text="Run MetCell",
                        command=run_metcell_subprocess, width=15)
btn_mc_run.pack(side='left', padx=5)

btn_mc_terminate = ttk.Button(btn_frame_mc, text="Terminate",
                              command=terminate_metcell_subprocess,
                              width=12, state="disabled")
btn_mc_terminate.pack(side='left', padx=5)

frame_metcell.rowconfigure(mc_row, weight=1)

# 初始化默认时间段
add_time_seg_row("0", "7200")
mc_vars["ion_mode"].trace('w', lambda *a: init_mc_bg_defaults())
init_mc_bg_defaults()

# ---------- 页：Cell-specific Peak Parameters（tims 模式）----------
cell_specific_tims_sf = ScrollableFrame(notebook)
frame_cs_tims = cell_specific_tims_sf.inner

row = 0
ttk.Label(frame_cs_tims, text="【Cell-specific peak filtration (IM-based)】",
          font=('Arial', 11, 'bold')).grid(
    row=row, column=0, columnspan=5, sticky='w', padx=5, pady=(10, 6))
row += 1

ttk.Label(frame_cs_tims, text="Use",       font=('Arial', 10, 'bold')).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs_tims, text="Parameter", font=('Arial', 10, 'bold')).grid(row=row, column=1, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="Threshold", font=('Arial', 10, 'bold')).grid(row=row, column=2, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="Lower",     font=('Arial', 10, 'bold')).grid(row=row, column=3, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="Default",   font=('Arial', 10, 'bold')).grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# SNR
mc_vars["use_snr"] = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs_tims, variable=mc_vars["use_snr"]).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs_tims, text="SNR:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["snr_thres"] = tk.StringVar(value="3.5")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["snr_thres"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
mc_vars["snr_lower"] = tk.StringVar(value="2")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["snr_lower"]).grid(row=row, column=3, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(3.5 / 2)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# Match ratio
mc_vars["use_match"] = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs_tims, variable=mc_vars["use_match"]).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs_tims, text="Match ratio:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["match_thres"] = tk.StringVar(value="0.67")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["match_thres"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
mc_vars["match_lower"] = tk.StringVar(value="0.3")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["match_lower"]).grid(row=row, column=3, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(0.67 / 0.3)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

# RSD subs
mc_vars["use_rsd"] = tk.BooleanVar(value=True)
ttk.Checkbutton(frame_cs_tims, variable=mc_vars["use_rsd"]).grid(row=row, column=0, padx=8, pady=3)
ttk.Label(frame_cs_tims, text="RSD subs:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["rsd_thres"] = tk.StringVar(value="0.2")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["rsd_thres"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
mc_vars["rsd_lower"] = tk.StringVar(value="0.1")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["rsd_lower"]).grid(row=row, column=3, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(0.2 / 0.1)", foreground="gray").grid(row=row, column=4, sticky='w', padx=8, pady=3)
row += 1

ttk.Separator(frame_cs_tims, orient='horizontal').grid(row=row, column=0, columnspan=5, sticky='ew', pady=10)
row += 1

# Score threshold
ttk.Label(frame_cs_tims, text="Score threshold:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["threshold"] = tk.StringVar(value="10")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["threshold"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(default 10)", foreground="gray").grid(row=row, column=3, sticky='w', padx=8, pady=3)
row += 1

# spec_num_lim
ttk.Label(frame_cs_tims, text="Limitation of detected spectra in cell pulse:").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["spec_num_lim"] = tk.StringVar(value="1")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["spec_num_lim"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(default 1, recommand 2 when cell pulse spectra number > 10)",
          foreground="gray").grid(row=row, column=3, sticky='w', padx=8, pady=3)
row += 1

# Cell occurrence
ttk.Label(frame_cs_tims, text="Cell occurrence (%):").grid(row=row, column=1, sticky='w', padx=8, pady=3)
mc_vars["cell_occurrence"] = tk.StringVar(value="5")
ttk.Entry(frame_cs_tims, width=10, textvariable=mc_vars["cell_occurrence"]).grid(row=row, column=2, sticky='w', padx=8, pady=3)
ttk.Label(frame_cs_tims, text="(default 5)", foreground="gray").grid(row=row, column=3, sticky='w', padx=8, pady=3)
row += 1

ttk.Separator(frame_cs_tims, orient='horizontal').grid(row=row, column=0, columnspan=5, sticky='ew', pady=10)
row += 1

# 按钮：Run MPCell + Terminate
btn_frame_cs_tims = ttk.Frame(frame_cs_tims)
btn_frame_cs_tims.grid(row=row, column=0, columnspan=5, pady=15)

btn_worker_tims_run = ttk.Button(btn_frame_cs_tims, text="Run MPCell",
                                 command=run_worker_tims_subprocess,
                                 width=15)
btn_worker_tims_run.pack(side='left', padx=10)

btn_cs_tims_terminate = ttk.Button(btn_frame_cs_tims, text="Terminate",
                                   command=terminate_worker_tims_subprocess,
                                   width=15, state="disabled")
btn_cs_tims_terminate.pack(side='left', padx=10)
row += 1

# ---------- 页 7：MetCell Output（tims 模式）----------
mc_output_sf = ScrollableFrame(notebook)
frame_mcoutput = mc_output_sf.inner

mc_result_text = tk.Text(frame_mcoutput, wrap=tk.WORD, height=60)
scrollbar_mc = ttk.Scrollbar(frame_mcoutput, command=mc_result_text.yview)
mc_result_text.configure(yscrollcommand=scrollbar_mc.set)
mc_result_text.pack(side='left', fill='both', expand=True, padx=5, pady=5)
scrollbar_mc.pack(side='right', fill='y')


# ============================================================
# 底部按钮
# ============================================================
button_frame = ttk.Frame(root)
button_frame.pack(side='bottom', fill='x', pady=10)

btn_quit = ttk.Button(button_frame, text="Exit", command=root.destroy, width=12)
btn_quit.pack(side='right', padx=5)


# ============================================================
# 页签重建 & 模式切换
# ============================================================
def update_bottom_buttons():
    """Run/Terminate 已移至 Cell-specific Peak Parameters 页，无需动态切换"""
    pass


def rebuild_notebook():
    for tab_id in notebook.tabs():
        notebook.forget(tab_id)

    notebook.add(file_settings_sf, text="1. File Settings")

    if mode_var.get() == "no_tims":
        notebook.add(param_setting_sf,  text="2. Parameter Setting")
        notebook.add(cell_specific_sf,  text="3. Cell-specific Peak Parameters")
        notebook.add(output_sf,         text="4. Output")
        notebook.add(quant_sf,          text="5. Quantification")
        notebook.add(anno_sf,           text="6. Metabolite Annotation")
    else:
        notebook.add(metcell_sf,   text="2. MetMPCell Parameter Setting")
        notebook.add(cell_specific_tims_sf, text="3. Cell-specific Peak Parameters")
        notebook.add(mc_output_sf, text="4. MetMPCell Output")

    notebook.select(0)


def on_mode_changed(*args):
    """模式切换：重建页签；no_tims 时同步刷新参数默认值"""
    if mode_var.get() == "no_tims":
        update_mode_params()
    rebuild_notebook()


# ============================================================
# 初始化流程
# ============================================================
ion_mode_var.trace('w', lambda *a: update_mode_params())
update_mode_params()
mode_var.trace('w', on_mode_changed)
rebuild_notebook()

root.mainloop()