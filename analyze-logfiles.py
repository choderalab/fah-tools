#!/usr/bin/env python

def parse_log(results_filename):
    """
    Parse the log file from a results packet.

    ARGUMENTS

    results_filename (string) - name of compressed results file to test

    RETURNS

    logtext - text of log file
    logdata - dict of important log contents

    """

    # Create temporary directory.
    import os, os.path, tempfile, shutil
    cwd = os.getcwd()
    tmpdir = tempfile.mkdtemp()

    # Extract source directory.
    [directory, filename] = os.path.split(results_filename)

    # Copy results to temporary directory.
    shutil.copyfile(results_filename, os.path.join(tmpdir, 'results.tar.bz2'))

    # Change to temporary directory.
    os.chdir(tmpdir)

    # Extract payload and results.
    import commands
    command = 'bzcat results.tar.bz2 | tar x' 
    commands.getoutput(command)
    
    # Read log file.
    log_filename = 'log.txt'
    logtext = read_file(log_filename)

    # Extract useful info from log file.
    logdata = dict()
    import re
    for line in logtext.split('\n'):
        m = re.match('^(.+?):(.+)', line)
        if m:
            groups = m.groups()
            key = groups[0].strip()
            value = groups[1].strip()
            logdata[key] = value
            # TODO: Add support for values that can span multiple lines, like Options and Args.

    # Clean up temporary directory.
    os.chdir(cwd)
    for filename in os.listdir(tmpdir):
        os.unlink(os.path.join(tmpdir, filename))
    os.removedirs(tmpdir)

    return (logtext, logdata)

example = """
**************************** Zeta Folding@home Core ****************************
       Type: 23
       Core: Zeta
    Website: http://folding.stanford.edu/
  Copyright: (c) 2009-2013 Stanford University
     Author: Yutong Zhao <yutong.zhao@stanford.edu>
       Args: -dir 03 -suffix 01 -version 703 -lifeline 20282 -checkpoint 15 -gpu
             1 -gpu-vendor nvidia
     Config: <none>
************************************ Build *************************************
    Version: 0.0.45
       Date: May 20 2013
       Time: 10:30:56
    SVN Rev: 4000
     Branch: fah/trunk/cores/gpu/fahcore
   Compiler: GNU 4.6.3
    Options: -std=gnu++98 -O3 -funroll-loops -mfpmath=sse -ffast-math
             -fno-unsafe-math-optimizations -msse2
   Platform: linux2 3.2.0-33-generic
       Bits: 64
       Mode: Release
************************************ System ************************************
        CPU: Intel(R) Xeon(R) CPU X5680 @ 3.33GHz
     CPU ID: GenuineIntel Family 6 Model 44 Stepping 2
       CPUs: 24
     Memory: 15.66GiB
Free Memory: 11.98GiB
    Threads: POSIX_THREADS
Has Battery: false
 On Battery: false
 UTC offset: -5
        PID: 20286
        CWD: /mnt/ramdisk/work
         OS: Linux 3.9.0-2-generic x86_64
    OS Arch: AMD64
       GPUs: 4
      GPU 0: NVIDIA:2 GF100 [GeForce GTX 480]
      GPU 1: NVIDIA:2 GF100 [GeForce GTX 480]
      GPU 2: NVIDIA:2 GF100 [GeForce GTX 480]
      GPU 3: NVIDIA:2 GF100 [GeForce GTX 480]
       CUDA: Not detected
********************************************************************************
[1] compatible platform(s):
  -- 0 --
  PROFILE = FULL_PROFILE
  VERSION = OpenCL 1.1 CUDA 4.2.1
  NAME = NVIDIA CUDA
  VENDOR = NVIDIA Corporation

(4) device(s) found on platform 0:
  -- 0 --
  DEVICE_NAME = GeForce GTX 480
  DEVICE_VENDOR = NVIDIA Corporation
  DEVICE_VERSION = OpenCL 1.1 CUDA

  -- 1 --
  DEVICE_NAME = GeForce GTX 480
  DEVICE_VENDOR = NVIDIA Corporation
  DEVICE_VERSION = OpenCL 1.1 CUDA

  -- 2 --
  DEVICE_NAME = GeForce GTX 480
  DEVICE_VENDOR = NVIDIA Corporation
  DEVICE_VERSION = OpenCL 1.1 CUDA

  -- 3 --
  DEVICE_NAME = GeForce GTX 480
  DEVICE_VENDOR = NVIDIA Corporation
  DEVICE_VERSION = OpenCL 1.1 CUDA

[ Entering Init ]
  Launch time: 2013.04.21  21:30:49
  Arguments passed: -dir 03 -suffix 01 -version 703 -lifeline 20282 -checkpoint 15 -gpu 1 -gpu-vendor nvidia 
[ Leaving  Init ]
[ Entering Main ]
  Reading core settings...
  Total number of steps: 2500000
  XTC write frequency: 50000
[ Initializing Core Contexts ]
  Using platform OpenCL
  Looking for vendor: nvidia...found on platformId 0
  Deserializing System...
  Setting up Force Groups:
    Group 0: Everything Else
    Group 1: Nonbonded Direct Space
    Group 2: Nonbonded Reciprocal Space
  Found MonteCarloBarostat @ 1.01325 (default) Bar, 300 Kelvin, 50 pressure change frequency.
    Found: 55480 atoms, 6 forces.
  Deserializing State...  done.
    Integrator Type: N6OpenMM18LangevinIntegratorE
    Constraint Tolerance: 1e-05
    Time Step in PS: 0.002
    Temperature: 300
    Friction Coeff: 5
  Checking core state against reference...
  Checking checkpoint state against reference...
[ Initialized Core Contexts... ]
  Using OpenCL on platformId 0 and gpu 1
  v(^_^)v  MD ready starting from step 0

...

2013.04.22  12:38:1
[ Leaving  Main ]
Saving result file logfile_01.txt
Saving result file checkpointState.xml
Saving result file checkpt.crc
Saving result file log.txt
"""

def histogram(logs, key):
    hist = dict()
    for log in logs:
        if key in log:
            value = log[key]
            if value in hist:
                hist[value] += 1
            else:
                hist[value] = 1
    return hist

def gpu_histogram(logs):
    keys = ['GPU 0', 'GPU 1', 'GPU 2', 'GPU 3']
    for key in keys:
        for log in logs:
            if key in log:
                value = log[key]
                elements = value.split()
                value = ' '.join(elements[1:])
                if value in hist:
                    hist[value] += 1
                else:
                    hist[value] = 1
    
    return hist

def generate_sunburst(filename, attribute, hist):
    outfile = open(filename, 'w')
    outfile.write('{\n')    
    outfile.write(' "name": "%s",\n' % attribute)
    outfile.write(' "children": [\n')
    for key in hist:
        value = hist[key]
        outfile.write('  {"name": "%s", "size": %d},\n' % (key, value))
    outfile.write(' ]\n')
    outfile.write('}\n')
    return
            
if __name__ == '__main__':
    # Read pickled log data.
    logdata_output_filename = '/cbio/jclab/projects/fah/output/logs.pkl'
    import cPickle
    logdata_outfile = open(logdata_output_filename, 'r')
    logs = cPickle.load(logdata_outfile)
    logdata_outfile.close()

    attributes = ['CPU', 'CPU ID', 'CPUs', 'Memory', 'Threads', 'Has Battery', 'On Battery']
    for attribute in attributes:
        hist = histogram(logs, attribute)

        print attribute
        print '-' * (48+1+8)
        for key in hist:
            print "%-48s %8d" % (key, hist[key])
        print ""
        
    # Pool GPU stats.
    hist = gpu_histogram(logs)

    print "GPUs (aggregated)"
    print '-' * (48+1+8)
    for key in hist:
        print "%-48s %8d" % (key, hist[key])
    print ""

    #import os.path
    #filename = os.path.join('output', attribute + '.html')
    #generate_sunburst(filename, hist)
        
    
