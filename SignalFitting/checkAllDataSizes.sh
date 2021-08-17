#!/bin/bash

#pt bins:
root -b -l -q "CheckDataSize.C(kAADATA, 0, 3, 0.0, 2.4, 10, 90, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 3, 6, 0.0, 2.4, 10, 90, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 6, 10, 0.0, 2.4, 10, 90, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 10, 50, 0.0, 2.4, 10, 90, 3.5)"

#y bins:
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 0.0, 0.8, 10, 90, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 0.8, 1.6, 10, 90, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 1.6, 2.4, 10, 90, 3.5)"

#centrality bins
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 0.0, 2.4, 10, 30, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 0.0, 2.4, 30, 50, 3.5)"
root -b -l -q "CheckDataSize.C(kAADATA, 0, 50, 0.0, 2.4, 50, 90, 3.5)"

