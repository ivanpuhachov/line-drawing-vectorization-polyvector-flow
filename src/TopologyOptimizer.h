#pragma once
#include "graph_typedefs.h"
#include "typedefs.h"
#include "FrameFieldFlow.h"
#include "ChainUtils.h"
#include "Coverage.h"
#include "SelectCurves.h"
#include "RemoveShortBranches.h"

std::array<std::pair<size_t, double>, 2> findClosestVerticesInReebGraph(const Eigen::Vector2d& p, int component, const G& reebGraph, const cv::Mat& origMask, const std::vector<int>& components, const std::array<Eigen::MatrixXcd, 2>& allRoots, const bool doTheGraphDistanceCheck);
std::vector<MyPolyline> optimizeTopology(G& reebGraph, const std::vector<cv::Point2d>& pts, const cv::Mat& bwImg, const cv::Mat& origMask, std::map<vertex_descriptor, bool>& covered, const FrameFieldFlow& fff, const std::array<Eigen::MatrixXcd, 2>& allRoots, std::vector<cv::Point2d>& newlyAddedIntersections, std::vector<std::vector<vertex_descriptor>>& finalchains, std::vector<std::vector<vertex_descriptor>>& allowedChains, std::map<int , int>& keypointToValenceMap, std::vector<std::pair<int, int>>& allowedChainToLabeledPointPair);
std::vector<cv::Point2d> addJunctions(G& g, const std::set<size_t>& activeVertices, const std::vector<cv::Point2d>& pts);