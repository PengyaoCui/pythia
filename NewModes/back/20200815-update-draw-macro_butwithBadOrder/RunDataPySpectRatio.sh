#!/bin/bash
#
for i in \
"\"Kshort\"" \
"\"Lambda_sum\"" \
"\"Xi\"" \
"\"Omega\"" \

do

echo .q | sleep 5 | root -l DataPySpectRatio.C\(${i}\)
done
#
exit 0

