#!/bin/env python

"""
Automated benchmarking of Folding@home projects

"""

import time
import datetime
import copy
import re
import numpy as np

projectid = 10490 # project number to benchmark
log_filename = "log.txt" # log file to process

def parse_logfile(filename, verbose=False):
    """
    Parse the specified FAHClient logfile.

    Parameters
    ----------
    filename : str
        The log filename to parse.
    verbose : bool, optional, default=False
        If True, print verbose output.

    Returns
    -------
    gpunames : list of str
        The list of GPU names for which benchmarking data is available.
    statistics : dict of dict
        statistics[gpuname] a dict of timing statistics (in days) for the given GPU.
        statistics[gpuname]['TPP'] is time per percent, in days
        statistics[gpuname]['TPPstd'] is standard deviation of the time per percent, in days
        statistics[gpuname]['TPWU'] is time per work unit, in days
        statistics[gpuname]['TPWUstd'] is standard deviation of the time per work unit, in days        
        
    """
    if verbose: print "Reading log file '%s'..." % filename

    # CONSTANTS
    SECONDS_PER_DAY = float(24*60*60) # number of seconds in a day

    # Read log file.
    infile = open(log_filename, 'r')
    lines = infile.readlines()
    infile.close()

    # Get reference date/time.
    # 
    # *********************** Log Started 2015-08-22T17:59:19Z ***********************
    result = re.match('^\*+ Log Started (\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\dZ) \*+\n', lines[0])
    logstart_datetime = result.groups()[0]
    logstart_struct_time = time.strptime(logstart_datetime, '%Y-%m-%dT%H:%M:%SZ')
    logstart_seconds = time.mktime(logstart_struct_time) # time since epoch in seconds at which log was started

    # Parse logfile.
    gpuname_in_slot = dict() # gpuname_in_slot[slot] is the name of the GPU in slot `slot`
    events = dict() # events[gpuname] is a list of all percent-completed reports, in seconds since start of log
    for line in lines[1:]:
        # Split into timestamp and message
        timestamp, message = line[0:8], line[9:]
        # Convert into struct_time.
        hour, min, sec = re.split(':',timestamp)
        # Compute elapsed time (in seconds) since log start.
        elapsed_seconds = datetime.timedelta(seconds=float(sec)-logstart_struct_time.tm_sec, 
                                             minutes=float(min)-logstart_struct_time.tm_min, 
                                             hours=float(hour)-logstart_struct_time.tm_hour).total_seconds()
        # Handle midnight rollover.
        if (elapsed_seconds < 0):
            elapsed_seconds += 24*60*60
    
        # Truncate message.
        message = message.strip()

        # Find GPU definitions
        # 17:59:19:Enabled folding slot 00: READY gpu:0:GK110 [GeForce GTX 780]
        result = re.match("^Enabled folding slot (\d+): READY gpu:\d:\S+ \[(.+)\]$", message)
        if result is not None:
            slot, gpuname = result.groups()
            if verbose: print "Found GPU '%s' in folding slot %s" % (gpuname, slot)
            gpuname_in_slot[slot] = gpuname

        # Find reported percentages
        # WU00:FS00:0x18:Completed 0 out of 5000000 steps (0%)
        result = re.match("^WU(\d+):FS(\d+):0x(\d+):Completed \d+ out of \d+ steps \((\d+)\%\)$", message)
        if result is not None:
            WU, slot, core, percent = result.groups()
            gpuname = gpuname_in_slot[slot]
            if gpuname in events:
                events[gpuname].append(elapsed_seconds)
            else:
                events[gpuname] = [elapsed_seconds]

    # Extract GPU names.
    slots = gpuname_in_slot.keys()
    gpunames = [ gpuname_in_slot[slot] for slot in slots ]

    if verbose:
        for gpuname in gpunames:
            print "%-20s: completed %d%%" % (gpuname, len(events[gpuname])-1)

    # Compute time per percent samples for each GPU
    percent_completed_timings = dict() # percent_completed_timings[gpuname] is a numpy array of all timing statistics for 1%
    for gpuname in gpunames:
        nevents = len(events[gpuname])
        percent_completed_timings[gpuname] = np.zeros([nevents-1], np.float64)
        for index in range(0,nevents-1):
            percent_completed_timings[gpuname][index] = events[gpuname][index+1]-events[gpuname][index]
        
    # Compute statistics.
    statistics = dict()
    for gpuname in gpunames:
        statistics[gpuname] = dict()
        statistics[gpuname]['TPP'] = percent_completed_timings[gpuname].mean() / SECONDS_PER_DAY # days
        statistics[gpuname]['TPPstd'] = percent_completed_timings[gpuname].std() / SECONDS_PER_DAY # days
        statistics[gpuname]['TPWU'] = statistics[gpuname]['TPP'] * 100
        statistics[gpuname]['TPWUstd'] = statistics[gpuname]['TPPstd'] * 100

    return gpunames, statistics

#
# REFERENCE BENCHMARK CARDS
#

benchmarks = { 
    'GeForce GTX 780' : 175000, # GTX-780 should get 175000 PPD
    }


def compute_benchmark(statistics, benchmarks=benchmarks, Kfactor=0.75):
    """
    Apply benchmarks using recommended PPD for a given set of cards.
    
    Parameters
    ----------
    statistics : dict of dict
        statistics[gpuname] a dict of timing statistics (in days) for the given GPU.
    benchmarks : dict of float
        benchmarks[gpuname] is the expected PPD for a given fixed benchmark card
    Kfactor : float, optional, default=0.75
        The bonus point factor

    Returns
    -------
    statistics : dict of dict
        statistics[gpuname] is the input statistics[gpuname] augmented with some additional statistics
        statistics[gpuname]['basepoints'] is the base points
        statistics[gpuname]['timeout'] is the timeout (in days)
        statistics[gpuname]['deadline'] is the deadline (in days)

    """

    for gpuname in statistics.keys():
        if gpuname in benchmarks:
            TPP = statistics[gpuname]['TPP'] # time per percent (in days)
            PPD = benchmarks[gpuname] # expected PPD
            timeout = 3.5 * statistics[gpuname]['TPWU'] # timeout (in days)
            deadline = timeout + 1 # deadline (in days)
            basecredit = float(PPD*1e3) / np.sqrt(Kfactor * deadline / TPP**3)
            
            statistics[gpuname]['timeout'] = timeout
            statistics[gpuname]['deadline'] = deadline
            statistics[gpuname]['basecredit'] = basecredit

    return statistics


if __name__ == '__main__':
    # Parse the log file.
    (gpunames, statistics) = parse_logfile(log_filename, verbose=True)

    # Compute benchmarks
    statistics = compute_benchmark(statistics)

    print ""
    print "BENCHMARK RESULTS"
    for gpuname in gpunames:
        if 'basecredit' in statistics[gpuname]:
            print "%-20s | %10.8f days/WU | timeout %8.6f deadline %8.6f credit %8f" % (gpuname, statistics[gpuname]['TPWU'], statistics[gpuname]['timeout'], statistics[gpuname]['deadline'], statistics[gpuname]['basecredit'])
        else:
            print "%-20s | %10.8f days/WU" % (gpuname, statistics[gpuname]['TPWU'])
    print ""


