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
* `set-pme-parameters-in-system.py` - script to add explicit PME parameters to system.xml files lacking them

## Usage

Set paths in `extract-logs.py`, `extract-logs-multiprocessing.py`, `analyze-logfiles.py`
```
python extract-logs-multiprocessing.py
python analyze-logfiles.py
```


## Benchmark Procedure (Updated 9/29/15)
* Login to `choderalab@csk.mskcc.org`
* `cd fah-7.4.9-centos`
* `rm -rf logs/ cores/ work/`
* Edit `config.xml` to have your project key 
* run `./FAHClient` for a few percent of the work unit 
* run `python ./fah-tools/benchmark.py` 
* Take points from GTX 780 
