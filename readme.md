# MPCellMet

**MPCellMet** is a multi-parameter-assisted, GUI-based data processing platform for
single-cell metabolomics (SCM) acquired by continuous-infusion mass cytometry.
At its core is the **MPCell framework**, which identifies cell-specific mass peaks
against background and noise signals by combining three complementary parameters —
signal-to-noise ratio (SNR), match ratio (MR) and RSD subs — through a piecewise
scoring framework with per-parameter rigid and loose thresholds.

MPCellMet supports two acquisition modes:

| Mode | Input | Upstream peak detection | Notes |
|---|---|---|---|
| **no-IM** | `.mzML` | built-in EIC extraction | for conventional continuous-infusion mass cytometry |
| **IM-based** | Bruker `.d` folder | MetCell (via Docker) | for ion-mobility-resolved mass cytometry |

Both modes share the same MPCell scoring framework, parameter semantics and output
format, so results are directly comparable across platforms.

---

## Contents

- [Features](#features)
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Parameters](#parameters)
- [Output files](#output-files)
- [Reproducing the reported performance](#reproducing-the-reported-performance)
- [Third-party components](#third-party-components)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

---

## Features

- **Explicit cell-specific peak identification.** Most existing SCM workflows retain
  all detected peaks or filter only by cell occurrence frequency. MPCellMet evaluates
  each candidate peak with three physically interpretable parameters and a piecewise
  score, which preserves rare, low-occurrence cell-specific features that
  frequency-only filtering discards.
- **Interpretable by design.** No model pre-training and no annotated training set are
  required. Every parameter can be adjusted per sample, per peak and per instrument,
  and every score is traceable to a measurable quantity.
- **Cross-platform compatibility.** Validated on Q Exactive Plus across both
  ionization polarities, nine cultured cell lines, primary mouse BMSCs and sorted
  human colorectal cancer cells, and simulated cell-pulse spans of 3–15 spectra.
- **Integrated GUI.** A Windows executable with task-oriented pages for
  parameter setting, peak filtering, quantification and metabolite annotation.
  No programming experience required.
- **Seamless MetCell integration.** For IM-based data, MPCellMet automatically
  launches the MetCell workflow inside a Docker container, then hands the candidate
  peaks to the MPCell framework — no R installation or manual format conversion.

---

## System requirements

| Component | Requirement |
|---|---|
| Operating system | Windows 10 / 11 (64-bit) |
| Python (source install only) | 3.9.x |
| MSConvert (no-IM mode) | ProteoWizard 3.0.x (tested with 3.0.22354-a648f68) |
| Docker (IM-based mode) | Docker Desktop 4.x with WSL2 backend |
| Disk | ≥ 5 GB free (Docker image + intermediate results) |
| RAM | ≥ 8 GB recommended |

The pre-built `MPCellMet.exe` bundles Python and all Python dependencies; users
running the executable do **not** need to install Python.

---

## Installation

### Option A — Pre-built Windows executable (recommended for end users)

1. Download the latest `MPCellMet-vX.Y.Z-win64.zip` from the
   [Releases](https://github.com/Zhangyi-pku/MPCell/releases) page.
2. Extract the archive to a folder of your choice.
3. (no-IM mode) Install [ProteoWizard MSConvert](http://proteowizard.sourceforge.net/)
   and make sure `msconvert.exe` is available (either on `PATH` or specified in the GUI).
4. (IM-based mode) Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
   and start it. Then pull the MetCell image once:
   ```bash
   docker pull zhulab/metcell-r:latest
   ```
   The image used in the published work is pinned to
   `zhulab/metcell-r@sha256: 493f11f12180cc386bf05335d51f53b514d29c5a60cbed2fad7a4c49e4842c3e`
   (pulled on 2026-11-20).

5. Double-click `MPCellMet.exe`.

### Option B — From source (for developers)

```bash
git clone https://github.com/Zhangyi-pku/MPCell.git
cd MPCell
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python MPCell_main_20261125.py
```

`requirements.txt` pins all Python dependencies:

```
numpy==1.24.4
pandas==2.0.3
matplotlib==3.7.5
pyteomics==4.6.3
openpyxl==3.1.2
```

`tkinter` is part of the Python standard library. If you build your own executable,
pin `pyinstaller==5.13.2` (the version used for the released `MPCellMet.exe`).

**MetCell (IM-based mode).** R and all MetCell dependencies are provided inside the
Docker image; users do **not** need to install R. `libtimsdata.so` is bundled with the
MetCell image and is used under Bruker's own license terms.

---

## Quick start

### no-IM mode

1. Convert your raw file to centroid mzML with MSConvert (MS1 only).
2. Launch `MPCellMet.exe` → select **no_IM** on the *File Settings* page.
3. On the *Parameter Setting* page, specify:
   - input `.mzML` file and output directory;
   - ion polarity;
   - peak intensity threshold, m/z tolerance and RT range;
   - cell marker m/z and intensity threshold;
   - one or more background markers (m/z + maximum intensity).
4. On the *Cell-specific Peak Parameters* page, keep the default thresholds
   (SNR 3.5/2, MR 0.67/0.3, RSD subs 0.2/0.1; score cutoff 10; cell occurrence 5%)
   or adjust them for your instrument.
5. Click **Run**. Results are written to the output directory; the *Output* page
   streams the log in real time.

### IM-based mode

1. Launch `MPCellMet.exe` → select **IM-based** on the *File Settings* page.
2. On the *MetMPCell Parameter Setting* page, specify:
   - Bruker `.d` folder and output directory;
   - one or more acquisition time segments (start/end, in seconds);
   - ion mode, marker m/z and mobility, m/z and mobility tolerances,
     intensity upper/lower limits, background markers.
3. Click **Run MetCell**. MPCellMet will:
   1. write `time_limit_table.csv`;
   2. create a Docker volume and copy the `.d` folder, MetCell R scripts and
      `libtimsdata.so` into it;
   3. run `Rscript /data/run.R` inside the MetCell image;
   4. copy results back to the output directory;
   5. automatically launch the MPCell module to extract EICs and identify
      cell-specific peaks.
4. The *MetMPCell Output* page shows progress and errors; a **Terminate** button is
   available on both the MetCell and MPCell pages.

---

## Parameters

Full definitions and default values are given in **Table S1** of the Supporting
Information. The most important parameters are:

| Parameter | Default | Meaning |
|---|---|---|
| Peak intensity threshold | `1e4` | minimum MS1 intensity for a peak to enter EIC construction |
| m/z tolerance | `10 ppm` | tolerance for EIC extraction, marker/background matching and clustering |
| Marker m/z (pos / neg) | `760.5851` / `306.0767` | cell marker used to detect cell pulses |
| Marker intensity threshold | `1e6` / `5e6` | minimum marker intensity defining a cell pulse |
| Cell EIC peak span | `5` | symmetric window around each marker apex |
| Limited cell EIC peak span | `1` | minimum number of consecutive above-threshold spectra for a cell pulse (≥2 recommended when a pulse has >10 spectra) |
| `n_spec` | `1` | minimum consecutive above-threshold spectra for a peak to count as present in a pulse |
| SNR rigid / loose | `3.5` / `2.0` | |
| MR rigid / loose | `0.67` / `0.30` | |
| RSD subs rigid / loose | `0.20` / `0.10` | |
| Cell occurrence threshold | `5%` | minimum fraction of cell pulses for a peak to be retained |
| Score threshold | `10` | cutoff on the summed score `S = S_SNR + S_MR + S_RSD` |

**Gating note.** MR uses the fraction of cell pulses with non-zero intensity (> 5%)
as an additional gate, while SNR and RSD subs use the cell occurrence ratio
(COR > 5%). The final decision always requires `S ≥ 10` **and** `COR ≥ 5%`.

---

## Output files

For a no-IM run, the output directory contains:

| Path | Content |
|---|---|
| `output.txt` | full processing log |
| `marker_intensity.xlsx` | marker EIC used for cell-pulse detection |
| `cell_marker.xlsx` | detected cell pulses (RT of first spectrum, max marker intensity) |
| `mz_tar.npy` | candidate peak m/z list after clustering |
| `c_peak_eic_intensity.npy` | EIC matrix (peaks × time) |
| `grouping/RT_del/cell_specific_peak.xlsx` | **cell-specific peaks** with SNR, MR, RSD subs, COR and score |
| `grouping/RT_del/cell_nonspecific_peak.xlsx` | non-specific peaks |
| `grouping/RT_del/all_candidate_peak.xlsx` | **all candidate peaks** with scores (useful for manual curation and benchmarking) |
| `grouping/RT_del/cell_specific_peak.npy` | array of cell-specific m/z values |
| `plot_cell-specific/` | EIC plots (optional, controlled by the *EIC Plot* setting) |

For an IM-based run, the same files are produced after MetCell completes, plus the
MetCell intermediate results in `results/`.

---

## Reproducing the reported performance

The example dataset below reproduces the FPR/FNR reported in the paper.

**Files** (available from the [Releases](https://github.com/Zhangyi-pku/MPCell/releases)
page and the Zenodo archive):

- `example/MDA-MB-231_pos_example.mzML` — representative no-IM acquisition.
- `example/example_params.json` — parameter file for this dataset.
- `example/example_reference_peak_list.xlsx` — manually annotated reference peaks,
  built from the **complete candidate peak list** (not from MPCell output) by two
  independent annotators; disagreements were resolved by a third annotator.
  Columns: `mz`, `manual_label` (`true`/`false`), `annotator_1`, `annotator_2`,
  `consensus`, `notes`.
- `example/example_MPCell_output.xlsx` — `all_candidate_peak.xlsx` produced by
  MPCellMet on the same file, for direct comparison.

**Steps**

1. Run MPCellMet on `MDA-MB-231_pos_example.mzML` with `example_params.json`.
2. Open `grouping/RT_del/all_candidate_peak.xlsx` and
   `example_reference_peak_list.xlsx`.
3. Merge the two tables on `mz` (within the m/z tolerance).
4. Compute:
   ```
   TP = peaks with manual_label = true  and score ≥ threshold (and COR ≥ 5%)
   FP = peaks with manual_label = false and score ≥ threshold (and COR ≥ 5%)
   FN = peaks with manual_label = true  and score <  threshold (or COR <  5%)
   FPR = FP / (FP + TP)
   FNR = FN / (FN + TP)
   ```
5. The resulting FPR and FNR should match those reported in the paper within
   rounding error.

**Annotated reference set independence.** The reference set was constructed from the
full candidate list, before any MPCell filtering. It is therefore independent of the
algorithm and can be used to evaluate sensitivity without circularity.

---

## Third-party components

MPCellMet depends on the following third-party packages and resources, each
distributed under its own license:

| Component | License |
|---|---|
| [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/) | BSD 3-Clause |
| [pyteomics](https://pyteomics.readthedocs.io/) | Apache License 2.0 |
| [openpyxl](https://openpyxl.readthedocs.io/) | MIT |
| tkinter / Tcl/Tk | Python Software Foundation License |
| [ProteoWizard MSConvert](http://proteowizard.sourceforge.net/) | Apache License 2.0 (used externally) |
| MetCell (`zhulab/metcell-r` Docker image) | see MetCell repository; used via Docker, not redistributed |
| Bruker `libtimsdata.so` | Bruker Daltonics license; bundled inside the MetCell image |

---

## License

MPCellMet is released under the **MIT License**. See [LICENSE](LICENSE) for the full
text.

```
Copyright (c) 2026 Yi Zhang, Mingying Shi, Mingxuan Li, Yu Bai
```

---

## Citation

If you use MPCellMet in your work, please cite:

> Zhang, Y.; Shi, M.; Li, M.; Bai, Y. *MPCellMet: A Multi-Parameter-Assisted
> Single-Cell Metabolomics Data Processing Platform.* (2026).
> DOI: `10.5281/zenodo.<XXXXXXX>`

A versioned, citable archive of the code corresponding to the published paper is
available on Zenodo at the DOI above. Please cite the specific version you used.

---

## Contact

- Corresponding author: **Yu Bai** — yu.bai@pku.edu.cn
- Issues and feature requests: please open a
  [GitHub issue](https://github.com/Zhangyi-pku/MPCell/issues).

## Needs update: 10.5281/zenodo.<XXXXXXX>（Zenodo 发布后的 DOI）、vX.Y.Z（你的发布版本号）