#include "graph_typedefs.h"
#include <boost/graph/filtered_graph.hpp>

template <typename Graph>
Graph makeSubGraph(
        const Graph& g,
        std::vector<vertex_descriptor> vertices
        ) {
    Graph subgraph(boost::num_vertices(g));
    for (size_t i=0; i<boost::num_vertices(g); ++i) {
        // copying vertex properties
        Cluster c = g[i];
//        c.location = Eigen::Vector2d(point.y, point.x);
//        c.split = false;
//        c.clusterIdx = reebGraph[vertexToMove].clusterIdx;
//        c.width = reebGraph[vertexToMove].width;
//        c.solvedWidth = reebGraph[vertexToMove].solvedWidth;
        subgraph[i] = c;
    }
    edge_iter eit, eend;
    for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
    {
        if (std::find(vertices.begin(), vertices.end(), eit->m_target) != vertices.end()) {
            auto e = boost::add_edge(eit->m_source, eit->m_target, subgraph);
            // copying edge properties
            subgraph[e.first].edgeCurve = g[*eit].edgeCurve;
            subgraph[e.first].weight = g[*eit].weight;
        }
    }
    return subgraph;
}

template <typename Graph>
std::vector<std::vector<typename Graph::edge_descriptor>> chainDecomposition(
        const Graph& g,
        std::map<typename Graph::edge_descriptor, size_t>& myChain
        )
{
	std::vector<std::vector<typename Graph::edge_descriptor>> chains;
	//1. Break into chains, record their adjacencies into another graph

	for (auto [eit, eend] = boost::edges(g); eit != eend; ++eit)
	{
		//find a seed edge for a new chain
		if (myChain.find(*eit) != myChain.end())
			continue;

		//grow the chain into both direction until it hits a high-valence vertex or reaches an end
		std::vector<typename Graph::edge_descriptor> newChain = { *eit };
		myChain[*eit] = chains.size();


		for (int dir : { -1, 1 })
		{
			size_t curVtx = dir == -1 ? newChain.front().m_source : newChain.back().m_target;
			bool continuationFound = true;
			//try extending starting from this vertex
			while (continuationFound && (boost::degree(curVtx, g) == 2))
			{
				continuationFound = false;
				
				for (auto [oeit, oeend] = boost::out_edges(curVtx, g); oeit != oeend; ++oeit)
				{
					if (myChain.find(*oeit) == myChain.end())
					{
						curVtx = oeit->m_target;
						if (dir == -1)
							newChain.insert(newChain.begin(), boost::edge(oeit->m_target, oeit->m_source, g).first);
						else
							newChain.push_back(*oeit);

						myChain[*oeit] = chains.size();
						continuationFound = true;
						break;
					}
				}
			}
		}

		chains.push_back(newChain);
	}
	return chains;
}

template <typename Graph>
std::vector<std::vector<typename Graph::edge_descriptor>> chainDecompositionWithActive(
        const Graph& g,
        std::set<size_t> activepoints,
        std::map<typename Graph::edge_descriptor, size_t>& myChain
)
{
    std::vector<std::vector<typename Graph::edge_descriptor>> chains;
    //1. Break into chains, record their adjacencies into another graph

    for (auto [eit, eend] = boost::edges(g); eit != eend; ++eit)
    {
        //find a seed edge for a new chain
        if (myChain.find(*eit) != myChain.end())
            continue;

        //grow the chain into both direction until it hits a high-valence vertex or reaches an end
        std::vector<typename Graph::edge_descriptor> newChain = { *eit };
        myChain[*eit] = chains.size();


        for (int dir : { -1, 1 })
        {
            size_t curVtx = dir == -1 ? newChain.front().m_source : newChain.back().m_target;
            bool continuationFound = true;
            //try extending starting from this vertex
            while (continuationFound && (boost::degree(curVtx, g) == 2) && (activepoints.find(curVtx) == activepoints.end()))
            {
                continuationFound = false;

                for (auto [oeit, oeend] = boost::out_edges(curVtx, g); oeit != oeend; ++oeit)
                {
                    if (myChain.find(*oeit) == myChain.end())
                    {
                        curVtx = oeit->m_target;
                        if (dir == -1)
                            newChain.insert(newChain.begin(), boost::edge(oeit->m_target, oeit->m_source, g).first);
                        else
                            newChain.push_back(*oeit);

                        myChain[*oeit] = chains.size();
                        continuationFound = true;
                        break;
                    }
                }
            }
        }

        chains.push_back(newChain);
    }
    return chains;
}

#ifdef __GNUG__
template <typename Graph, typename EdgePredicate, typename VertexPredicate>
std::vector<std::vector<typename boost::filtered_graph<Graph,EdgePredicate,VertexPredicate>::edge_descriptor>> chainDecomposition(const boost::filtered_graph<Graph,EdgePredicate,VertexPredicate>& g, std::map<typename boost::filtered_graph<Graph,EdgePredicate,VertexPredicate>::edge_descriptor, size_t>& myChain)
{
    // hack to make GCC use valid chainDecomposition for filtered graph
    return chainDecomposition<boost::filtered_graph<Graph,EdgePredicate,VertexPredicate>>(g, myChain);
}
#endif