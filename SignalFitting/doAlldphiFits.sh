#!/bin/bash
whichSyst=0
whichRound=0

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

#integrated bin:
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
#pt bins:
./DoFit.sh kAADATA 0.0 3.0 0.0 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 3.0 6.0 0.0 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 6.0 10.0 0.0 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 10.0 50.0 0.0 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
#y bins:
./DoFit.sh kAADATA 0.0 50.0 0.0 0.8 10 90 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.8 1.6 10 90 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 1.6 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 1.2 10 90 3.5 0.00 0.12 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 1.2 2.4 10 90 3.5 0.00 0.12 $whichSyst $whichRound
#centrality bins
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 30 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 50 3.5 0.00 0.12 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 50 90 3.5 0.00 0.12 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 50 3.5 0.00 0.12 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 90 3.5 0.00 0.12 $whichSyst $whichRound

#integrated bin:
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
#pt bins:
./DoFit.sh kAADATA 0.0 3.0 0.0 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 3.0 6.0 0.0 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 6.0 10.0 0.0 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 10.0 50.0 0.0 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
#y bins:
./DoFit.sh kAADATA 0.0 50.0 0.0 0.8 10 90 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.8 1.6 10 90 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 1.6 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 1.2 10 90 3.5 0.12 0.25 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 1.2 2.4 10 90 3.5 0.12 0.25 $whichSyst $whichRound
#centrality bins
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 30 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 50 3.5 0.12 0.25 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 50 90 3.5 0.12 0.25 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 50 3.5 0.12 0.25 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 90 3.5 0.12 0.25 $whichSyst $whichRound

#integrated bin:
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
#pt bins:
./DoFit.sh kAADATA 0.0 3.0 0.0 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 3.0 6.0 0.0 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 6.0 10.0 0.0 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 10.0 50.0 0.0 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
#y bins:
./DoFit.sh kAADATA 0.0 50.0 0.0 0.8 10 90 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.8 1.6 10 90 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 1.6 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 1.2 10 90 3.5 0.25 0.38 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 1.2 2.4 10 90 3.5 0.25 0.38 $whichSyst $whichRound
#centrality bins
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 30 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 50 3.5 0.25 0.38 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 50 90 3.5 0.25 0.38 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 50 3.5 0.25 0.38 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 90 3.5 0.25 0.38 $whichSyst $whichRound

#integrated bin:
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
#pt bins:
./DoFit.sh kAADATA 0.0 3.0 0.0 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 3.0 6.0 0.0 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 6.0 10.0 0.0 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 10.0 50.0 0.0 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
#y bins:
./DoFit.sh kAADATA 0.0 50.0 0.0 0.8 10 90 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.8 1.6 10 90 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 1.6 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 1.2 10 90 3.5 0.38 0.5 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 1.2 2.4 10 90 3.5 0.38 0.5 $whichSyst $whichRound
#centrality bins
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 30 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 50 3.5 0.38 0.5 $whichSyst $whichRound
./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 50 90 3.5 0.38 0.5 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 50 3.5 0.38 0.5 $whichSyst $whichRound
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 30 90 3.5 0.38 0.5 $whichSyst $whichRound

