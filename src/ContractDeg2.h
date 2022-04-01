#pragma once
#include "graph_typedefs.h"

template <typename Graph>
std::vector<std::pair<typename Graph::vertex_descriptor, typename Graph::vertex_descriptor>> topoGraph(const Graph& g)
{
	std::map<edge_descriptor, size_t> myChain;
	auto chains = chainDecomposition(g, myChain);
	std::vector<std::pair<typename Graph::vertex_descriptor, typename Graph::vertex_descriptor>> result(chains.size());
	for (int i = 0; i < chains.size(); ++i)
		result[i] = { chains[i].front().m_source,chains[i].back().m_target };
	return result;
}

std::vector<std::pair<size_t, size_t>> topoGraphHighValenceSeparated(const G & g, std::vector<std::vector<edge_descriptor>>& chainsSeparated, bool onlyLoops=false);
void contractDeg2(G& g);
void contractAllButChainEnds(G & g, const std::vector<std::vector<vertex_descriptor>>& finalchains);
