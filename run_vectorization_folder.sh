#!/bin/sh

# This script runs vectorization for each .png and .jpg file in examples/

for pngfile in examples/*.png
do
    ptsfile=${pngfile//.png/.pts}
    svgfile=${pngfile//.png/_vector.svg}
    txtfile=${pngfile//.png/_vector_logs.txt}

	echo "\n\n\n ${pngfile} - $ptsfile - $svgfile - $txtfile"
    > ${txtfile}
    bash predict_points.sh $pngfile $ptsfile | tee ${txtfile}

    cmake-build-release/vectorize ${pngfile} $ptsfile $svgfile | tee ${txtfile}
done



for pngfile in examples/*.jpg
do
    ptsfile=${pngfile//.jpg/.pts}
    svgfile=${pngfile//.jpg/_vector.svg}
    txtfile=${pngfile//.jpg/_vector_logs.txt}

	echo "\n\n\n ${pngfile} - $ptsfile - $svgfile - $txtfile"
    > ${txtfile}
    bash predict_points.sh $pngfile $ptsfile | tee ${txtfile}
    
    cmake-build-release/vectorize ${pngfile} $ptsfile $svgfile | tee ${txtfile}
done
