#!/bin/bash

mkdir -p $3/$4

# blender 3.3
$blender=blender

# results dir
$results=/cluster/home/threnggli/results


for i in $(seq -f "%06g" $1 $2); do
    echo rendering $i
    blender -b automatic.blend --python automatic.py -o $3/$4/img_######.png -f $i -- $results/$3/interpolated/mesh-$i.ply > /dev/null
done

cd $3/$4 
ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
