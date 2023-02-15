#!/usr/bin/env bash

seeds=(11 3459 1346653 246 38479)
make
sed -i "/seed/s/= .*/= ${seeds[0]}/" input.toml
./main
for i in {1..4}
do
    mv energy_data.out energy_data_$i.out
    mv position_data.out position_data_$i.out
    sed -i "/seed/s/= .*/= ${seeds[$i]}/" input.toml
    ./main
done
./src/analysis.py
