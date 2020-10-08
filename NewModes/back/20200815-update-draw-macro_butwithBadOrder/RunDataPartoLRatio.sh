#!/bin/bash
#
for i in \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l dataPartoLRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l dataPartoLRatio_pPb.C\(${i},1\)
echo .q | sleep 5 | root -l dataPartoLRatio_pPb.C\(${i},0\)
done
#
exit 0

