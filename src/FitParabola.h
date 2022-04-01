#ifndef _FIT_PARABOLA_H_
#define _FIT_PARABOLA_H_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "opencv2/core/core.hpp"

bool fitParabola(const cv::Mat & bwImg, float strokeRadius, const cv::Mat & mask, double x0, double y0, bool onlyReliable, Eigen::Vector2d& outCenterPt, Eigen::Vector2d& outCenterAxis, std::array<double, 2> & outD, double & outResidual, int & outNotNull);

#endif
