#!/bin/bash
#
for i in \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 2 | root -l draw_DataPy.C\(${i}\)
done
#
exit 0

