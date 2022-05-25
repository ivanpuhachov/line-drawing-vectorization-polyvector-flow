#!/bin/sh

# This script runs vectorization pipeline on imagefile provided (.png or .jpg)

image=${1:-"examples/dog06.png"}

ptsfile=${image//.png/_auto.pts}
svgfile=${image//.png/_filter.svg}
txtfile=${image//.png/_logs.txt}

ptsfile=${ptsfile//.jpg/_auto.pts}
svgfile=${svgfile//.jpg/_filter.svg}
txtfile=${txtfile//.jpg/_logs.txt}

echo $image - $ptsfile - $svgfile

bash predict_points.sh $pngfile $ptsfile

cmake-build-release/polyvector_thing $image $ptsfile $svgfile | tee $txtfile