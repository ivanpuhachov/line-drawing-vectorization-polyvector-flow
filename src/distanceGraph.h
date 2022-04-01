#include "graph_typedefs.h"
#include "typedefs.h"

/***
 * This function computes graph edge distances for each point in polys to interesting polys
 * @param reebGraph
 * @param polys
 * @param allowedChains
 * @return
 */
std::vector<std::vector<double>> computeGraphDistancesPerPoly(const G& reebGraph, std::vector<MyPolyline> polys, const std::vector<std::vector<vertex_descriptor>>& allowedChains);

