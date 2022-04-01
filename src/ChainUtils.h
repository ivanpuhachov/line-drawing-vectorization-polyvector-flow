#pragma once
#include "typedefs.h"
#include "graph_typedefs.h"
#include "intersections.h"

void printChain(std::vector<edge_descriptor>c, const G& g, bool onlyEnds=false, bool printCoords=false);

void printChainByEdges(std::vector<edge_descriptor>c, G& g);

MyPolyline chainToPolyline (const Chain& c, const G& g);

MyPolyline vertexSeqToPolyline(const std::vector<vertex_descriptor>& s, const G& g);

double chainLength(std::vector<edge_descriptor> c, const G& g);

template <typename Graph> std::vector<typename Graph::edge_descriptor> flipChain(const Graph& g, const std::vector<typename Graph::edge_descriptor>& chain);

Chain mergeChains (Chain chainI, Chain chainJ, const G& g);

bool areChainsIntersecting(std::vector<edge_descriptor> c1, std::vector<edge_descriptor> c2, G& g);

std::vector<vertex_descriptor> chainToVerticesSeq(const Chain& c);

bool vertexIsInsideChain(const vertex_descriptor& v, const Chain& chain);

std::pair<Chain, Chain> splitChainByVertex(const vertex_descriptor& v, const Chain& chain);

int getVertexValence(const vertex_descriptor& v, const std::vector<Chain>& chains);

Chain alignChainToStartWith(const G& g, const vertex_descriptor& v, Chain& chain);

bool chainsOverlapFromStart(const Chain& chain1, const Chain& chain2);

std::pair<std::pair<Chain, Chain>, std::pair<Chain, Chain>> decomposeOverlappingChains(const Chain& chain1, const Chain& chain2);

bool checkIfEdgeIsInChains (const edge_descriptor& edge, const std::vector<std::vector<edge_descriptor>>& chains);

vertex_descriptor chainElementStartingFrom(
        const std::vector<edge_descriptor>& selectedChain,
        const size_t startElement,
        const size_t idx
);

vertex_descriptor chainElementStartingFrom(
        const std::vector<std::vector<edge_descriptor>>& chains,
        const int chainIdx,
        const size_t startElement,
        const size_t idx
);

/**
 * Returns pair of chain ends
 * @param chains
 * @param chainIdx
 * @return
 */
std::pair<vertex_descriptor, vertex_descriptor> getChainEndsIndicies (
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx
);
/**
 * Returns the other end of a chain
 * @param chains
 * @param chainIdx
 * @param knownEnd
 * @return
 */
vertex_descriptor getOtherEndOfChain (
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx,
        size_t knownEnd);

/**
 * This function estimates the coverage of a chain as a sum of solvedWidths radiuses
 * @param g
 * @param chains
 * @param chainIdx
 * @param startingVertexIdx
 * @return
 */
double estimateCoverageOfAChainStartingFrom(
        const G& g,
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx,
        size_t startingVertexIdx);

/**
 * This function splits all of the chains such that no active point is in the middle of a chain.
 * @param chain
 * @param activePts
 * @return
 */
std::vector<std::vector<edge_descriptor>>
splitChainAtActivePts(
        std::vector<edge_descriptor>chain,
        std::set<std::pair<int, int>> activePts);