#!/bin/bash
#
for i in \
"\"Kshort\"","kTRUE" \
"\"Kshort\"","kFALSE" \
"\"Lambda_sum\"","kTRUE" \
"\"Lambda_sum\"","kFALSE" \
"\"Xi\"","kTRUE" \
"\"Xi\"","kFALSE" \
"\"Omega\"","kTRUE" \
"\"Omega\"","kFALSE" \

do

echo .q | sleep 5 | root -l PtSpect_DataPy.C\(${i}\)
done
#
exit 0

