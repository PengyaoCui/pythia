#!/bin/bash
#
for i in \
"\"Xi\"","kTRUE" \
"\"Omega\"","kTRUE" \
"\"Xi\"","kFALSE" \
"\"Omega\"","kFALSE" \

do

echo .q | sleep 5 | root -l DataPyPartoLRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l DataPyPartoLRatio_pPb.C\(${i}\)
done
#
exit 0

