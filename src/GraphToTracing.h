#ifndef _GRAPH_TO_TRACING_H_
#define _GRAPH_TO_TRACING_H_
#include "typedefs.h"
#include "graph_typedefs.h"

void graphToTracing(G& g, const cv::Mat& origMask, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<std::vector<edge_descriptor>>& chains, bool secondSolve);

#endif