#pragma once

#include "typedefs.h"
#include <set>
double chopThisCurve(const std::vector<MyPolyline>& polys, int i, int whichEnd, const std::vector<std::array<bool, 2>>& protectedEnds, double threshold);