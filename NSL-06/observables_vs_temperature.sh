#!/usr/bin/env bash

./clean.sh

while read temp; do
    echo "### Computing $temp T"
    sed -i .bak "1s/TEMPERATURE/$temp/" input.dat
    ./Monte_Carlo_ISING_1D.exe
    sed -i .bak "1s/$temp/TEMPERATURE/" input.dat
done <temperatures.dat

# Nota: per utilizzo su GNU togliere .bak
