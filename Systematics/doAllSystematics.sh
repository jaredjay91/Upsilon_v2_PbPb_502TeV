#!/bin/bash

root -l -b -q "getSystematics.C(1,1)"
root -l -b -q "getSystematics.C(1,2)"
root -l -b -q "getSystematics.C(1,3)"
root -l -b -q "getSystematics.C(1,4)"

root -l -b -q "mergeSystematics.C()"
