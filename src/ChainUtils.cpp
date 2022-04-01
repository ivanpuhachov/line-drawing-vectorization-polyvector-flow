#include "ChainUtils.h"


void printChain(std::vector<edge_descriptor>c, const G& g, bool onlyEnds, bool printCoords){
    if (onlyEnds) {
        std::cout<<" ( " <<c.front().m_source<< " -> " << c.back().m_target << " ) " << std::endl;
    }
    else
    {
        std::cout << "edge start [ ";
        for (auto e : c)
        {
            std::cout<<e.m_source;
            if (printCoords)
                std::cout << " ( "<<g[e.m_source].location.x()<<", "<<g[e.m_source].location.y() << " )";
            std::cout<<" - ";
        }
        std::cout<<c.back().m_target;
        if (printCoords)
            std::cout<< " ( "<<g[c.back().m_target].location.x()<<", "<<g[c.back().m_target].location.y()<<" )";
        std::cout<<" ] edge end "<<std::endl;
    }
};

void printChainByEdges(std::vector<edge_descriptor>c, G& g){
    for (auto e : c)
    {
        std::cout<<g[e.m_source].location.x()<<", "<<g[e.m_source].location.y()<<"  to  "<<g[e.m_target].location.x()<<", "<<g[e.m_target].location.y()<<std::endl;
    }
    std::cout<<"index : "<<c.back().m_target<< " Coords : "<<g[c.back().m_target].location.x()<<", "<<g[c.back().m_target].location.y()<<std::endl<<std::endl;
};

MyPolyline chainToPolyline (const Chain& c, const G& g)
{
    MyPolyline poly;
    for (const G::edge_descriptor& e : c)
    {
        poly.push_back(g[e.m_source].location);
    }
    poly.push_back(g[c.back().m_target].location);
    return poly;
};

MyPolyline vertexSeqToPolyline(const std::vector<vertex_descriptor>& s, const G& g){
    MyPolyline poly;
    for (const vertex_descriptor & v : s)
    {
        poly.push_back(g[v].location);
    }
    return poly;
}

double chainLength(std::vector<edge_descriptor> c, const G& g)
{
    double l = 0;
    for (auto e : c)
    {
        l += (g[e.m_target].location - g[e.m_source].location).norm();
    }
    return l;
}

template <typename Graph>
std::vector<typename Graph::edge_descriptor> flipChain(const Graph& g, const std::vector<typename Graph::edge_descriptor>& chain)
{
    std::vector<typename Graph::edge_descriptor> result;
    for (const auto& e : chain)
        result.push_back(boost::edge(e.m_target, e.m_source, g).first);

    std::reverse(result.begin(), result.end());
    return result;
}

Chain mergeChains (Chain chainI, Chain chainJ, const G& g)
{
    if (chainI.back().m_target != chainJ.front().m_source)
    {
        //need to reverse things a bit
        if (chainI.back().m_target == chainJ.back().m_target)
        {
            chainJ = flipChain(g, chainJ);
        }
        else if (chainI.front().m_source == chainJ.back().m_target)
        {
            std::swap(chainI, chainJ);
        }
        else if (chainI.front().m_source == chainJ.front().m_source)
        {
            chainI = flipChain(g, chainI);
        }
    }

    std::vector<edge_descriptor> mergedChain = chainI;
    mergedChain.insert(mergedChain.end(), chainJ.begin(), chainJ.end());
    return mergedChain;
};

bool areChainsIntersecting(std::vector<edge_descriptor> c1, std::vector<edge_descriptor> c2, G& g)
{
    for (auto& e1 : c1)
    {
        for (auto& e2 : c2)
        {
            double s, t; // WN : if intersection s returns ratio of first segment to intersection point, t same for 2nd
            auto p0 = g[e1.m_source].location, p1 = g[e1.m_target].location;
            auto p2 = g[e2.m_source].location, p3 = g[e2.m_target].location;

            if (get_line_intersection(p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), nullptr, nullptr, &s, &t))
                if ((s > 1e-6) && (s < 1 - 1e-6) && (t > 1e-6) && (t < 1 - 1e-6)) // WN : I guess that means that edges intersect ?
                    return true;
        }
    }
    return false;
}

std::vector<vertex_descriptor> chainToVerticesSeq(const Chain& c)
{
    std::vector<vertex_descriptor> result;
    for (const G::edge_descriptor& e : c)
    {
        result.push_back(e.m_source);
    }
    result.push_back(c.back().m_target);
    return result;
}

bool vertexIsInsideChain(const vertex_descriptor& v, const Chain& chain){
    /*
     * Returns true if vertex is INSIDE the chain (but not at the ends)
     */
    bool result = false;
    for (size_t i=0; i<chain.size()-1; ++i){
        if (chain[i].m_target==v){
            result = true;
        }
    }
    return result;
}

std::pair<Chain, Chain> splitChainByVertex(const vertex_descriptor& v, const Chain& chain){
    Chain chain1, chain2;
    bool buildingChain1 = true;
    for (auto edge : chain) {
        if (buildingChain1)
            chain1.push_back(edge);
        else
            chain2.push_back(edge);
        if (edge.m_target==v)
            buildingChain1 = false;
    }
    return std::make_pair(chain1, chain2);
}

int getVertexValence(const vertex_descriptor& v, const std::vector<Chain>& chains){
    int result = 0;
    for (auto & chain : chains){
        result += (chain.front().m_source==v);
        result += (chain.back().m_target==v);
    }
    return result;
}

Chain alignChainToStartWith(const G& g, const vertex_descriptor& v, Chain& chain){
    if (chain[0].m_source==v){
        return chain;
    } else {
        assert(chain.back().m_target==v);
        return flipChain(g, chain);
    }
}

bool chainsOverlapFromStart(const Chain& chain1, const Chain& chain2) {
    if (chain1.front().m_source!=chain2.front().m_source){
        std::cout << "ERROR: chain overlapping FAILS\n";
    }
    return (chain1.front().m_target==chain2.front().m_target);
}

std::pair<std::pair<Chain, Chain>, std::pair<Chain, Chain>> decomposeOverlappingChains(const Chain& chain1, const Chain& chain2) {
    size_t minlen = std::min(chain1.size(), chain2.size());
    size_t separatorId = 0;
    for (size_t edge_i=0; edge_i<minlen; ++edge_i){
        if (chain1[edge_i].m_target!=chain2[edge_i].m_target){
            separatorId=edge_i;
            break;
        }
    }
    assert(separatorId!=0);
    Chain overlap(chain1.begin(), chain1.begin()+separatorId);
    Chain residual1(chain1.begin()+separatorId, chain1.end());
    Chain residual2(chain2.begin()+separatorId, chain2.end());
    return std::make_pair(
            std::make_pair(overlap, residual1),
            std::make_pair(overlap, residual2)
    );
}

bool checkIfEdgeIsInChains (const edge_descriptor& edge, const std::vector<std::vector<edge_descriptor>>& chains){
    bool result = false;
    for (auto chain : chains) {
        result = result || (std::find(chain.begin(), chain.end(), edge) != chain.end());
    }
    return result;
};

vertex_descriptor chainElementStartingFrom(
        const std::vector<edge_descriptor>& selectedChain,
        const size_t startElement,
        const size_t idx
)
{
    if (selectedChain.front().m_source == startElement) // WN : The first point of the first edge of the chain == startElement ?
    {
        if (selectedChain.size() == idx) // WN : Fct is only called with idx = 1 or 2. Checks if the chain has 1 or 2 edges.
            return selectedChain.back().m_target; // WN : Then returns the last vertex

        return selectedChain[idx].m_source; // WN : Then return the second vertex of the chain (if idx = 1) or the third vertex of the chain (if idx = 2)
    }
    else // WN : The first vertex of the chain is not the startElement. So it is the last. Now if idx = 1 return the second to last vertex. If idx = 2 return the third to last vertex.
    {
        if (selectedChain.size() == idx) // WN : Fct is only called with idx = 1 or 2. So this checks if the chain is only 1 or 2 edges (depending on value of idx)
            return selectedChain.front().m_source; // WN : Then return the first vertex of the first edge of the chain.

        return selectedChain[selectedChain.size() - 1 - idx].m_target; // WN : Else return before last vertex (if idx = 1) or the before before last vertex (if idx = 2)
    }
}

vertex_descriptor chainElementStartingFrom(
        const std::vector<std::vector<edge_descriptor>>& chains,
        const int chainIdx,
        const size_t startElement,
        const size_t idx
)
{
    return chainElementStartingFrom(chains[chainIdx], startElement, idx);
};


std::pair<vertex_descriptor, vertex_descriptor> getChainEndsIndicies (
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx
){
    return std::make_pair(chains[chainIdx].front().m_source, chains[chainIdx].back().m_target);
};


vertex_descriptor getOtherEndOfChain (
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx,
        size_t knownEnd)
{
    std::pair ends = getChainEndsIndicies(chains, chainIdx);
    if (knownEnd == ends.first)
        return ends.second;
    if (knownEnd == ends.second)
        return ends.first;
    // if passed index is not on of the ends of specifiend chain throw exception
    throw std::invalid_argument("You tried to find the other end of a chain you don't belong to. Don't do this!");
};


double estimateCoverageOfAChainStartingFrom(
        const G& g,
        const std::vector<std::vector<edge_descriptor>>& chains,
        int chainIdx,
        size_t startingVertexIdx)
{
    // estimating the coverage of a chain WITHOUT the coverage of starting vertex
    std::vector<edge_descriptor> selectedChain = chains[chainIdx];
    double chainCoverage = - g[startingVertexIdx].solvedWidth;
    for (auto edge : selectedChain) {
        chainCoverage += g[edge.m_source].solvedWidth;
    }
    chainCoverage += g[selectedChain.back().m_target].solvedWidth;
    return chainCoverage;
};

std::vector<std::vector<edge_descriptor>> splitChainAtActivePts(std::vector<edge_descriptor>chain, std::set<std::pair<int, int>> activePts)
{
    std::vector<std::vector<edge_descriptor>> chains;

    std::vector<edge_descriptor> tmp_c;
    for(auto& e : chain)
    {
        bool foundAHit = false;
        if (e != chain.back())
        {
            for (auto &actPt : activePts) {
                if (e.m_target == actPt.second)
                {
                    std::cout << "found a hit ! ! " << std::endl;
                    foundAHit = true;
                    break;
                }
            }
        }
        tmp_c.push_back(e);
        if (foundAHit)
        {
            chains.push_back(tmp_c);
            tmp_c.clear();
            foundAHit = false;
        }
        if (e == chain.back())
            chains.push_back(tmp_c);
    }
    return chains;
}