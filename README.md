# MPCell

Windows 64-bit executable package.

## Download
Download `MPCell_v1.0.0_win64.zip` from [Releases](../../releases).
No need to download Souce code(zip) and Source code(tar.zpip).

## Usage
1. Extract the zip to a folder with an English-only path, e.g. `D:\MPCell\`.
2. Keep all `.exe` files and the `metcell_files/` folder together.
3. Double-click the desired `MPCell_main.exe` to run.

> Do not move the `.exe` files out of the folder, or they may fail to find dependencies.

Third-party components
----------------------
This software depends on the following third-party packages, each distributed
under its own license:

- pyteomics (Apache License 2.0)
- numpy, pandas, matplotlib (BSD 3-Clause)
- openpyxl (MIT)
- tkinter / Tcl/Tk (Python Software Foundation License)
- ProteoWizard MSConvert (Apache License 2.0), used externally for raw-to-mzML conversion
- MetCell (https://github.com/...) is used via its Docker image `zhulab/metcell-r`
  and is distributed under its own license; MPCellMet does not redistribute
  MetCell source code.
- Bruker libtimsdata.so is provided by Bruker Daltonics and is subject to
  Bruker's own license terms; users are responsible for obtaining it legally.
