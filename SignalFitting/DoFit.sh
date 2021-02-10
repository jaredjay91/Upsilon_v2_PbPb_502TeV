#!/bin/bash
collId=kAADATA
ptLow=0.0
ptHigh=30.0
yLow=0.0
yHigh=2.4
cLow=10
cHigh=60
muPtCut=4.0
dphiEp2Low=0.0
dphiEp2High=0.5
whichSyst=0
whichRound=0
nTries=1

if [ "$1" != "" ]
then
  collId=$1
  shift
fi
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
  muPtCut=$1
  shift
fi
if [ "$1" != "" ]
then
  dphiEp2Low=$1
  shift
fi
if [ "$1" != "" ]
then
  dphiEp2High=$1
  shift
fi
if [ "$1" != "" ]
then
  whichSyst=$1
  shift
fi
if [ "$1" != "" ]
then
  whichRound=$1
  shift
fi
if [ "$1" != "" ]
then
  nTries=$1
  shift
fi

if [ $nTries -lt 2 ]
then
  echo root -b -q -l "CheckFitExists.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,$whichRound)"
  root -b -q -l "CheckFitExists.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,$whichRound)" >> junk
  isIt=$(tail -c 2 junk)
  rm junk
else
  isIt=1
fi

if [ $isIt -gt 0 ]
then
  if [ $nTries -lt 2 ]
  then
    echo "THE FIT EXISTS ALREADY! :)"
  fi
  echo root -b -q -l "CheckFitQuality.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,$whichRound)"
  root -b -q -l "CheckFitQuality.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,$whichRound)" >> junk
  ToF=$(tail -c 2 junk)
  rm junk
else
  echo "THE FIT DOES NOT EXIST YET! :O"
  echo
  ToF=0
fi

if [ $ToF -gt 0 ]
then
  echo "THE FIT PASSED THE QUALITY CHECK! :)"
  echo
else
  echo "THE FIT FAILED THE QUALITY CHECK! :("
  echo
  if [ $nTries -lt 6 ]
  then
    echo "STARTING TRY #$nTries"
    if [ $nTries -lt 2 ]
    then
      echo root -b -q -l "FitDataWithRandomSeeds.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,kTRUE,$whichRound)"
      root -b -q -l "FitDataWithRandomSeeds.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,kTRUE,$whichRound)"
    else
      echo root -b -q -l "FitDataWithRandomSeeds.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,kTRUE,$whichRound)"
      root -b -q -l "FitDataWithRandomSeeds.C($collId,$ptLow,$ptHigh,$yLow,$yHigh,$cLow,$cHigh,$muPtCut,$dphiEp2Low,$dphiEp2High,$whichSyst,kTRUE,$whichRound)"
    fi
    nTries=$(($nTries+1))
    ./DoFit.sh $collId $ptLow $ptHigh $yLow $yHigh $cLow $cHigh $muPtCut $dphiEp2Low $dphiEp2High $whichSyst $whichRound $nTries
  else
    nTries=$(($nTries-1))
    echo "GIVING UP AFTER $nTries TRIES :("
    root -b -q -l "WriteToLog.C()"
    echo
  fi
fi

