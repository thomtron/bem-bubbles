#!/bin/bash

cd build
dir=$1
#mkdir $dir

for ply in "$dir"/*.ply
do
    ./remeshing-for-spherical-harmonics-pinned $ply
done

python ../python_utils/spherical-harmonics/spherical-harmonics-pinned.py $dir/sh-coeffs.csv $dir/mesh*.csv

    