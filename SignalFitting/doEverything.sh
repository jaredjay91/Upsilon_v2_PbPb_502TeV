#!/bin/bash

./doAllFullPhiFits.sh 0 R4a #Nominal fits
./doAllFullPhiFits.sh 0 R4b #Alternate constraint fits
./doAllFullPhiFits.sh 1 R4a #AltSig fits
./doAllFullPhiFits.sh 2 R4a #AltBkg fits
./doAllFullPhiFits.sh 3 R4a #AltAcc fits
./doAllFullPhiFits.sh 4 R4a #AltEff fits
./doAlldphiFits.sh 0 R4a #Nominal dphi fits
./doAlldphiFits.sh 0 R4b #Alternate constraint dphi fits
./doAlldphiFits.sh 1 R4a #AltSig dphi fits
./doAlldphiFits.sh 2 R4a #AltBkg dphi fits
./doAlldphiFits.sh 3 R4a #AltAcc dphi fits
./doAlldphiFits.sh 4 R4a #AltEff dphi fits
