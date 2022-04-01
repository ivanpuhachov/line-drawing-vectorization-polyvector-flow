#pragma once
#include "graph_typedefs.h"
#include "typedefs.h"
#include <queue>
#include "FrameFieldFlow.h"
#include "Coverage.h"

bool yJunctionTest(const std::vector<vertex_descriptor>& curPath, vertex_descriptor newV, const G& g);

std::vector<edge_descriptor> findShortestPath_AtLeastNNotCoveredVerticesMNewPixels(const G& g, const std::vector<vertex_descriptor>& myComponent, const std::set<vertex_descriptor>& terminals, const std::map<vertex_descriptor,bool>& covered, const int N, const FrameFieldFlow& fff, const std::set<std::pair<int, int>> & existingCoverage, const int M, const cv::Mat& mask);
size_t countSharedCurves(vertex_descriptor vertex1, vertex_descriptor vertex2, const G& g);