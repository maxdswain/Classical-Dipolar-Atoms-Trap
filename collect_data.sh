#!/usr/bin/env bash

seeds=(11 3459 1346653 246 38479)
make
for i in {1..5}
do
    ./main
    mv energy_data.txt energy_data_$i.txt
    sed -i "/seed/s/= .*/= ${seeds[$i-1]}/" input.toml
done
./src/analysis.py
