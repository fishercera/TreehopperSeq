#!/usr/bin/bash

echo "Library" > col1_tmp
grep "Input Read Pairs" *.err | cut -d : -f 1 | sed "s/\.polya.*//g" >> col1_tmp

echo "Input Read Pairs" > col2_tmp
grep "Input Read Pairs" *.err | cut -d : -f 3 | sed "s/Both Surviving//g" >> col2_tmp

echo "Both Surviving" > col3_tmp
grep "Input Read Pairs" *.err | cut -d : -f 4 | sed "s/Forward Only Surviving//g" >> col3_tmp

echo "Forward Only Surviving" > col4_tmp
grep "Input Read Pairs" *.err | cut -d : -f 5 | sed "s/Reverse Only Surviving//g" >> col4_tmp

echo "Reverse Only Surviving" > col5_tmp
grep "Input Read Pairs" *.err | cut -d : -f 6 | sed "s/Dropped//g" >> col5_tmp

echo "Dropped" > col6_tmp
grep "Input Read Pairs" *.err | cut -d : -f 7 >> col6_tmp

paste col1_tmp col2_tmp col3_tmp col4_tmp col5_tmp col6_tmp > TrimmingSummary.txt
rm -rf *tmp
