#!/bin/bash
#
for i in \
"\"Kshort\"" \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l dataSpect_pp.C\(${i}\)
echo .q | sleep 5 | root -l dataSpect_pPb.C\(${i},0\)
echo .q | sleep 5 | root -l dataSpect_pPb.C\(${i},1\)
done
#
exit 0

