#pragma once
#include "typedefs.h"
#include "graph_typedefs.h"
#include "ChainUtils.h"
#include "Coverage.h"
#include "FrameFieldFlow.h"
#include <gurobi_c++.h>

std::vector<Chain> selectBestCurves(const G& g, const std::vector<Chain>& currentChains, const cv::Mat& origMask, const FrameFieldFlow& fff);
std::vector<Chain> selectAllCurves(const G& g, const std::vector<Chain>& currentChains, const cv::Mat& origMask, const FrameFieldFlow& fff);

std::vector<MyPolyline> selectBestPolylines(const G& g,
                                            const std::vector<MyPolyline>& currentPolys,
                                            const std::vector<MyPolyline>& initialPolys,
                                            const std::vector<std::vector<vertex_descriptor>>& chains,
                                            const FrameFieldFlow& fff,
                                            std::vector<std::pair<int, int>> allowedChainToLabeledPointPair,
                                            const std::vector<int>& pts_types,
                                            std::map<std::pair<int,int>, std::vector<size_t>> pixelToChainCoveringMap,
                                            std::vector<bool>& selection);