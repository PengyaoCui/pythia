#!/bin/bash
#
for i in \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l DataPyPartoKRatioRatio.C\(${i}\)
done
#
exit 0

