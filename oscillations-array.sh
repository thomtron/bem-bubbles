#!/bin/bash

radius=($(seq 60e-6 5e-6 65e-6))
pressure=($(seq 13e4 1e4 14e4))

cd build
parentdir=oscillation-results
mkdir $parentdir

for pre in "${pressure[@]}"
do
    for rad in "${radius[@]}"
    do
        dir=$parentdir/results-$pre-Pa-$rad-m
        mkdir $dir
        ./oscillations $rad $pre $dir

        for ply in "$dir"/*.ply
        do
            ./remeshing-for-spherical-harmonics $ply
        done

        ls $dir

        python ../python_utils/spherical-harmonics/spherical-harmonics.py $dir/sh-coeffs.csv $dir/mesh*.csv

    done
done
    