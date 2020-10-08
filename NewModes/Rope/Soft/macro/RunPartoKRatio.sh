#!/bin/bash
#
for i in \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 2 | root -l PartoKRatio.C\(${i}\)
done
#
echo .q | sleep 2 | root -l PartoKRatio_JE.C

exit 0

