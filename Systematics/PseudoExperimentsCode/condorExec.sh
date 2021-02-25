#!/bin/bash
collId=kAADATA
ptLow=0.0
ptHigh=50.0
yLow=0.0
yHigh=2.4
cLow=10
cHigh=90
whichSyst=1
iTrial=1

if [ "$1" != "" ]
then
  ptLow=$1
  shift
fi
if [ "$1" != "" ]
then
  ptHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  yLow=$1
  shift
fi
if [ "$1" != "" ]
then
  yHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  cLow=$1
  shift
fi
if [ "$1" != "" ]
then
  cHigh=$1
  shift
fi
if [ "$1" != "" ]
then
  whichSyst=$1
  shift
fi
if [ "$1" != "" ]
then
  iTrial=$1
  shift
fi
cd /afs/cern.ch/work/j/jjay/public/Upsilonv2/Systematics/PseudoExperiments
. /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh
echo root -b -q "FitPseudoDataSimultaneously.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$whichSyst,$iTrial)";
root -b -q "FitPseudoDataSimultaneously.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$whichSyst,$iTrial)";
