#! /bin/bash

export HOME=/uscms/home/nmirman
cd /uscms/home/nmirman/CMSSW_5_3_5/src
eval `scramv1 runtime -sh`

export WORKING_DIR=$_CONDOR_SCRATCH_DIR

cd /uscms/home/nmirman/metsig/
./EvalSig -n 0.01
