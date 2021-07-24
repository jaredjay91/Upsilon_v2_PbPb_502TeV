#!/bin/bash

root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(1,1)"
root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(1,2)"
root -l -b -q "getSystematics.C(1,3)"
root -l -b -q "getSystematics.C(1,4)"
root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(1,5)"

root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(2,1)"
root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(2,2)"
root -l -b -q "getSystematics.C(2,3)"
root -l -b -q "getSystematics.C(2,4)"
root -l -b -q "PseudoExperimentsCode/getSystematicsFromPseudoExps.C(2,5)"

root -l -b -q "mergeSystematics.C(1)"
root -l -b -q "mergeSystematics.C(2)"

