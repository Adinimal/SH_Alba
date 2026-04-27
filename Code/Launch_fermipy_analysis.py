#!/usr/bin/env python

"""
Launcher for SH_Alba fermipy analysis.
Submits Fermi-LAT spectral hardening jobs for any source listed in an input text file.

Usage:
    python Launch_QPO_gating_mod3.py J0522_fullmission_windows.txt

If no argument is given, it launches the default input file.
"""

import sys, os

# ---- CONFIGURATION ----
BASE_DIR = '/sdf/home/a/adinesh/SH_Alba/Code/'
ANALYSIS_SCRIPT = '%s/SH_fermipy_analysis_V1.py' % BASE_DIR
LOG_DIR = '/sdf/home/a/adinesh/fermi-user/logs/SH_Alba'

# Default input files if no arguments given

DEFAULT_FILES = [
    '%s/J0522_fullmission_windows.txt' % BASE_DIR,
]

# ------------------------

# Parse command-line arguments
if len(sys.argv) > 1:
    input_files = sys.argv[1:]
else:
    input_files = DEFAULT_FILES
    
    print("No input files specified. Launching the default input file.")

total_submitted = 0

for INPUT_FILE in input_files:
    # Resolve path: if relative, prepend BASE_DIR
    if not os.path.isabs(INPUT_FILE):
        INPUT_FILE = os.path.join(BASE_DIR, INPUT_FILE)

    if not os.path.exists(INPUT_FILE):
        print("ERROR: Input file not found: %s" % INPUT_FILE)
        continue

    print("\n" + "=" * 70)
    print("Processing: %s" % INPUT_FILE)
    print("=" * 70)

    myf = open(INPUT_FILE, 'r')
    mylines = myf.readlines()
    myf.close()

    n_submitted = 0
    for i in range(len(mylines)):
        if mylines[i].startswith('#') or mylines[i].strip() == '':
            continue

        myline = mylines[i].split()
        srcname   = myline[0]          # e.g. 4FGLJ1555.7+1111 or 4FGLJ2158.8-3013
        timestart = float(myline[1])   # MJD start
        timeend   = float(myline[2])   # MJD end
        ID        = myline[3]          # e.g. T-6_p05, P+3_p10
        

        # ---- Split source name for the analysis script ----
        # name1 = "4FGL", name2 = "J1555.7+1111" or "J2158.8-3013"
        name1 = srcname[0:4]
        name2 = srcname[4:]

        
        # ---- Build the command ----
        cmd = "%s %s %s %.2f %.2f %s" % (ANALYSIS_SCRIPT, name1, name2,
                                          timestart, timeend, ID)

        # Log files organized by source/type
        log_subdir = '%s/%s' % (LOG_DIR, srcname)
        os.makedirs(log_subdir, exist_ok=True)
        logfile = '%s/logfile_%s_%s.txt' % (log_subdir, srcname, ID)


        # ---- SLURM settings ----
        queue = '--time=120:00:00'
        hosts = '--ntasks=1'
        cores = '--cpus-per-task=1'
        memor = '--mem-per-cpu=8g'
        parti = '--partition=milano'
        acc   = '--account=fermi:users'

        cmd = 'sbatch %s %s %s %s %s %s -o %s %s' % (queue, hosts, cores, memor,
                                                        acc, parti, logfile, cmd)
        print("Running:: %s" % cmd)
        os.system(cmd)
        n_submitted += 1

    print("\nSubmitted %d jobs for %s" % (n_submitted, os.path.basename(INPUT_FILE)))
    total_submitted += n_submitted

print("\n" + "=" * 70)
print("TOTAL SUBMITTED: %d jobs" % total_submitted)
print("=" * 70)
