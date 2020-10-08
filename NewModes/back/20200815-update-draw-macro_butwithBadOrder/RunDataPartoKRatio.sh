#!/bin/bash
#
for i in \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l dataPartoKRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l dataPartoKRatio_pPb.C\(${i},1\)
echo .q | sleep 5 | root -l dataPartoKRatio_pPb.C\(${i},0\)
done
#
exit 0

