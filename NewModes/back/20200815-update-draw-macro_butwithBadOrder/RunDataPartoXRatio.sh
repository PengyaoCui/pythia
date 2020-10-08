#!/bin/bash
#
for i in \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l dataPartoXRatio_pp.C\(${i}\)
echo .q | sleep 5 | root -l dataPartoXRatio_pPb.C\(${i},1\)
echo .q | sleep 5 | root -l dataPartoXRatio_pPb.C\(${i},0\)
done
#
exit 0

