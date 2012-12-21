#! /bin/bash

export HOME=/uscms/home/nmirman
cd /uscms/home/nmirman/CMSSW_4_4_4/src
eval `scramv1 runtime -sh`

export WORKING_DIR=$_CONDOR_SCRATCH_DIR

cd /uscms/home/nmirman/metsig/
./DoFit -j $1
