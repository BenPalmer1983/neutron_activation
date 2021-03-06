#!/bin/bash

# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=1
export PROC_COUNT=4
export PWSCF_SCRATCH=/opt/scratch
export PWSCF_PP=/opt/pp
export PWSCF_CACHE=/opt/pwscf_cache
export PWSCF_BIN=/opt/qe/bin/pw.x


thisdir=$(pwd)
cd ../../
topdir=$(pwd)
./package.sh
cd $thisdir
python3 $topdir/neutrons.py input.in

