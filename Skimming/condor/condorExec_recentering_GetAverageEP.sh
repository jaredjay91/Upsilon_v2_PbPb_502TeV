#!/bin/bash
cd /afs/cern.ch/work/j/jjay/public/Upsilonv2/condor
root -b -l -q "SkimRDTree_recentering_GetAverageEP.C($1,$2)" 
