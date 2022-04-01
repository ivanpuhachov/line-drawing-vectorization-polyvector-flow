#include "stdafx.h"
#include "OrientRoots.h"
#include <boost/graph/breadth_first_search.hpp>

class bfs_visitor : public boost::default_bfs_visitor
{
public:
	bfs_visitor(std::vector<bool>& vis) :visited(vis) {}

	void discover_vertex(size_t v, const G& g) const
	{
		visited[v] = true;
	}

	template <typename E, typename Graph>
	void examine_edge(E& e, Graph& g)
	{
		if (g[e.m_source].root.dot(g[e.m_target].root) < 0)
		{
			g[e.m_target].root = -g[e.m_target].root;
		}
	}
private:
	std::vector<bool>& visited;
};

void orientRoots(G & g, const std::vector<std::vector<edge_descriptor>>& chains, const std::map<edge_descriptor,size_t>& myChain)
{
//BFS
	std::vector<bool> visited(boost::num_vertices(g),false);
	bfs_visitor vis(visited);
	bool foundRoot = true;
	while (foundRoot)
	{
		foundRoot = false;
		for (size_t i = 0; i < visited.size(); ++i)
		{
			if (!visited[i])
			{
				foundRoot = true;
				boost::breadth_first_search(g, i, boost::visitor(vis));
			}
		}
	}

	std::map<edge_descriptor, Eigen::Vector2d> edgeRoots;
//fine. now roots are consistently oriented. for each chain decide how to orient them properly according to the chain direction
	for (int i = 0; i < chains.size(); ++i)
	{
		int sumDot = 0;
		double sign = 1.0; //tracks if we need to flip root sign since chain contains a previously deg-1 vertex
		for (auto e : chains[i])
		{
			size_t ms = e.m_source, mt = e.m_target;
			Eigen::Vector2d root0 = g[ms].root;
			Eigen::Vector2d root1 = g[mt].root;

			if (root0.norm() > 100) //for the newly added vertices at the intersections of two strokes
				root0 = root1;
			if (root1.norm() > 100)
				root1 = root0;

			//first, compute the average of two roots
			Eigen::Vector2d avgRoot = (root0 + root1)*0.5;

			bool orientationOK = root1.dot(root0) > 0;
			if (!orientationOK)
			{
				//this can happen only at the beginning of the chain
				std::cout << "THIS IS NOT IMPLEMENTED" << std::endl;
				//assert(false);
				avgRoot = (root0 - root1)*0.5;
			}

			g[e].avgRoot = sign*avgRoot;

			//now see if we should orient everything forward or backward
			double d = g[e].avgRoot.dot(g[mt].location - g[ms].location);
			if (d > 0)
				sumDot++;
			else
				sumDot--;

			if (g[mt].sharpCorner)
				sign = -sign;
		}

		if (sumDot < 0)
		{
			for (auto e : chains[i])
				g[e].avgRoot = -g[e].avgRoot;
		}
	}
}
