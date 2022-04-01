#pragma once
#include "typedefs.h"
#include "graph_typedefs.h"

std::vector<std::vector<vertex_descriptor>> findAllPaths(const G& g, vertex_descriptor s, vertex_descriptor t);