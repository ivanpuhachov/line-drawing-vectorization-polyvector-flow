#!/bin/sh

# This script runs vectorization pipeline on imagefile provided (.png or .jpg)

image=${1:-"examples/dog06.png"}

ptsfile=${image//.png/_auto.pts}
svgfile=${image//.png/_result.svg}
txtfile=${image//.png/_logs.txt}

ptsfile=${ptsfile//.jpg/_auto.pts}
svgfile=${svgfile//.jpg/_result.svg}
txtfile=${txtfile//.jpg/_logs.txt}

echo $image - $ptsfile - $svgfile

bash predict_points.sh $image $ptsfile

build/vectorize $image $ptsfile $svgfile
