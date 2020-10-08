#!/bin/bash
#
for i in \
"\"Lambda_sum\"","kTRUE" \
"\"Xi\"","kTRUE" \
"\"Omega\"","kTRUE" \
"\"Lambda_sum\"","kFALSE" \
"\"Xi\"","kFALSE" \
"\"Omega\"","kFALSE" \

do

echo .q | sleep 5 | root -l PartoKRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l PartoKRatio_pPb.C\(${i}\)
done
#
exit 0

