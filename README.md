# fah-tools
Tools for analyzing Folding@home logfiles

## Manifest
* `extract-logs.py` - script to pickle logfiles and perform various QC checks
* `extract-logs-multiprocessing.py` - multiprocessing version of pickling of logfiles
* `analyze-logfiles.py` - analyze pickled logfiles for statistics
* `quality_control.py` - quality control checks of SiegeTank data
* `optimizepme.py` - PME optimization script from Peter Eastman
* `openMM_testscript.py` - check energies across platforms
* `benchmark.py` - draft benchmarking script

## Usage

Set paths in `extract-logs.py`, `extract-logs-multiprocessing.py`, `analyze-logfiles.py`
```
python extract-logs-multiprocessing.py
python analyze-logfiles.py
```
