#!/bin/sh

# This scripts runs point prediction for all images in 'examples/' subfolder

for pngfile in examples/*.png
do
    ptsfile=${pngfile//.png/.pts}
    bash predict_points.sh $pngfile $ptsfile
done

for pngfile in examples/*.jpg
do
    echo $pngfile
    ptsfile=${pngfile//.jpg/_auto.pts}
    bash predict_points.sh $pngfile $ptsfile
done
