#!/bin/bash

nevt=-1
todaysDate=20201105
flattenBinByBin=kTRUE

#root -l -q "SkimRDTree_raw_GetAverageQ.C($nevt,$todaysDate)"
#root -l -q "SkimRDTree_recentering_GetAverageEP.C($nevt,$todaysDate)"
root -l -q "SkimRDTree_flatten_GetResCor.C($nevt,$todaysDate,$flattenBinByBin)"
