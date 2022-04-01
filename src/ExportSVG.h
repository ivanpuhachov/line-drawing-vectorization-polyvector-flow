#pragma once
#include "typedefs.h"

//saves polylines to the "filename.svg". If specified, will load a raster image from bgFilename and add it to the background of the svg
void exportSVG(std::string filename, const std::vector<MyPolyline>& polys, double width, double height, std::string bgFilename = "");
