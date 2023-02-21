#!/usr/bin/env bash

SEEDS=(11 3459 1346653 246 38479)
make
for i in {1..4}
do
    sed -i "/seed/s/= .*/= ${SEEDS[$i]}/" input.toml
    ./main & 
done
wait
./src/analysis.py
