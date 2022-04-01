#pragma once

#include "typedefs.h"
#include "graph_typedefs.h"

std::set<std::pair<int, int>> computeCoverage(const G& g, const std::vector<vertex_descriptor>& path, const cv::Mat& mask, bool includeEnds=true, bool circleAround=true);
std::set<std::pair<int, int>> computeTotalCoverage(const G& g, const std::vector<std::vector<vertex_descriptor>>& existingPaths,const cv::Mat& mask);
std::set<std::pair<int, int>> computeTotalCoverageByVertices(const G& g, const std::map<vertex_descriptor,bool>& covered);
int computeAddedCoverage(const G& g, const std::set<std::pair<int, int>>& totalCoverage, const std::vector<vertex_descriptor>& newPath, const cv::Mat& mask);
void addCoverage(const G& g, const std::vector<vertex_descriptor>& newPath, std::set<std::pair<int, int>>& totalCoverage, const cv::Mat& mask);
std::map<std::pair<int,int>, std::vector<size_t>> buildPixelToChainsCoveringMap(const G& g, const std::vector<std::vector<vertex_descriptor>>& existingPaths, const cv::Mat& mask);