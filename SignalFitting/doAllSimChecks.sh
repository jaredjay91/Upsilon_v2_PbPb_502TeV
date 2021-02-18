#!/bin/bash
whichSyst=0
whichRound=0

if [ "$1" != "" ]
then
  whichSyst=$1
  shift
fi

#integrated bin:
#./DoFit.sh kAADATA 0.0 50.0 0.0 2.4 10 90 3.5 $whichSyst 
#pt bins:
./DoSimCheck.sh kAADATA 0.0 3.0 0.0 2.4 10 90 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 3.0 6.0 0.0 2.4 10 90 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 6.0 10.0 0.0 2.4 10 90 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 10.0 50.0 0.0 2.4 10 90 3.5 $whichSyst 
#y bins:
./DoSimCheck.sh kAADATA 0.0 50.0 0.0 0.8 10 90 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 0.0 50.0 0.8 1.6 10 90 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 0.0 50.0 1.6 2.4 10 90 3.5 $whichSyst 
#centrality bins
./DoSimCheck.sh kAADATA 0.0 50.0 0.0 2.4 10 30 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 0.0 50.0 0.0 2.4 30 50 3.5 $whichSyst 
./DoSimCheck.sh kAADATA 0.0 50.0 0.0 2.4 50 90 3.5 $whichSyst 

