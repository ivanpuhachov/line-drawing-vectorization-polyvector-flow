#pragma once

#include "typedefs.h"

typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> Segment;
double perpendicularDistance(const Eigen::Vector2d& point, Segment segment, double& t);
MyPolyline simplify(const MyPolyline& poly, double eps);
MyPolyline simplify_lin_upsample(const MyPolyline& poly, double eps);
MyPolyline simplify_lin_downsample(const MyPolyline& poly, double eps);

