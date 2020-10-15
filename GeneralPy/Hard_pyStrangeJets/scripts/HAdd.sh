#!/bin/bash

PtHad=(005 011 021 036 057 "INF")
ptHad=(005 011 021 036 057)

if [ ! -d hard ]; then
  mkdir hard 
fi

#

for i in ${!ptHad[@]}
 
do
     hadd hard/AnalysisResults_PtHat_${PtHad[i]}_${PtHad[i+1]}.root ./dPtHat_${PtHad[i]}_${PtHad[i+1]}/AnalysisResults_*.root

done
#
exit 0

