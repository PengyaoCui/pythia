#!/bin/bash
#
for i in \
"\"Omega\"","kTRUE" \
"\"Omega\"","kFALSE" \

do

echo .q | sleep 5 | root -l DataPyPartoXRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l DataPyPartoXRatio_pPb.C\(${i}\)
done
#
exit 0

