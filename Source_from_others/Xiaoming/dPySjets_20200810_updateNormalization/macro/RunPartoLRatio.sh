#!/bin/bash
#
for i in \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l PartoLRatio.C\(${i}\)
done
#
exit 0

