#!/bin/sh

# This script runs vectorization for each .png and .jpg file in examples/

for pngfile in examples/*.png
do
    ptsfile=${pngfile//.png/.pts}
    svgfile=${pngfile//.png/_test_auto.svg}
    txtfile=${pngfile//.png/_auto_logs.txt}
    echo "${pngfile}- $ptsfile - $svgfile"
    cmake-build-release/vectorize ${pngfile} $ptsfile $svgfile | tee ${txtfile}
done



for pngfile in */*.jpg
do
    ptsfile=${pngfile//.jpg/_auto.pts}
    svgfile=${pngfile//.jpg/_auto.svg}
    txtfile=${pngfile//.jpg/_auto_logs.txt}
	  echo "${pngfile}- $ptsfile - $svgfile"
    cmake-build-release/vectorize ${pngfile} $ptsfile $svgfile | tee ${txtfile}
done
