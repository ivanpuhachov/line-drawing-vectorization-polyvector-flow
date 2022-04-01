#include "stdafx.h"
#include "GraphToTracing.h"
#include "Equations.h"
//#include <Eigen/SparseQR>
//#include <Eigen/IterativeLinearSolvers>
#include "GurobiLSLinIneqSolver.h"
#include <numeric>
#include <opencv2/highgui/highgui.hpp>
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "Params.h"

void graphToTracing(G & g, const cv::Mat& origMask, const std::array<Eigen::MatrixXcd, 2>& roots, const std::vector<std::vector<edge_descriptor>>& chains, bool secondSolve)
{
	std::cout << "Final optimization...";

	int n = boost::num_vertices(g);
	Equations eqs;

	//hold the vertices close to where they are
	const double inertiaWeight = OPTIMIZATION_INERTIA_WEIGHT;

	std::vector<double> weight(n, 1.0);
	std::map<edge_descriptor, double> dist;
	{
		edge_iter eit, eend;
		for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
			dist[*eit] = (g[eit->m_target].location - g[eit->m_source].location).norm();
	}

	for (size_t v = 0; v < n; ++v)
	{
		if (g[v].somethingHappened)
		{
			double r = g[v].width;
			
			std::vector<vertex_descriptor> pDij(n);
			std::vector<double> dDij(n);
			auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, g));
			auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx, g));
			auto fixedDistMap = boost::const_associative_property_map< std::map<edge_descriptor, double> >(dist);
			dijkstra_shortest_paths(g, v,
				predecessor_map(predMap).
				distance_map(distMap).weight_map(fixedDistMap));

			for (size_t v1 = 0; v1 < n; ++v1)
			{
				double d = distMap[v1];
				double sigma = 2 * r;
				double w = 1-exp(-d*d / (2 * sigma*sigma));
				weight[v1] = std::min(weight[v1], w);
			}
		}
	}

	for (size_t i = 0; i < n; ++i)
	{
		bool ignore = false;
		boost::graph_traits<G>::out_edge_iterator eit, eend;
	
		for (std::tie(eit, eend) = boost::out_edges(i, g); eit != eend; ++eit)
		{
			if (g[*eit].weight < 1e-5)
				ignore = true;
		}

		Eigen::Vector2d p0 = g[i].location;
		std::array<int, 2> p0i = { (int)std::round(p0.y()), (int)std::round(p0.x()) };

		if (ignore || (p0i[0] >= origMask.rows) || (p0i[0] < 0) || (p0i[1] >= origMask.cols) || (p0i[1] < 0) || (origMask.at<uchar>(p0i[0], p0i[1]) == 0))
			continue;

		double w = weight[i];//1.0;

		Eigen::Vector2d perp(-g[i].root.y(), g[i].root.x());

		eqs.addEqLHS(i * 2, inertiaWeight*w);
		eqs.addEqRHSAndFinishEquation(inertiaWeight*p0.x()*w);

		eqs.addEqLHS(i * 2 + 1, inertiaWeight*w);
		eqs.addEqRHSAndFinishEquation(inertiaWeight*p0.y()*w);
		//eqs.addIneqLHS(2 * i, inertiaWeight*w*perp.x());
		//eqs.addIneqLHS(2 * i + 1, inertiaWeight*w*perp.y());
		//eqs.addIneqRHSAndFinishEquation(inertiaWeight*w*perp.dot(g[i].location));

		//eqs.addIneqLHS(2 * i, inertiaWeight*w*g[i].root.x());
		//eqs.addIneqLHS(2 * i + 1, inertiaWeight*w*g[i].root.y());
		//eqs.addIneqRHSAndFinishEquation(inertiaWeight*w*g[i].root.dot(g[i].location));

	}

	//parallelism
	boost::graph_traits<G>::edge_iterator eit, eend;
	size_t eIdx = 0;
	double parallelismWeight = OPTIMIZATION_PARALLELISM_WEIGHT;
	for (int i = 0; i < chains.size(); ++i)
	{
		for (auto e: chains[i])
		{
			size_t ms = e.m_source, mt = e.m_target;
			
			double w = 1;
			eqs.addEqLHS(mt * 2, parallelismWeight*w);
			eqs.addEqLHS(ms * 2, -parallelismWeight*w);
			eqs.addEqLHS(2 * n + eIdx, -parallelismWeight*g[e].avgRoot.x()*w);
			eqs.addEqRHSAndFinishEquation(0);

			eqs.addEqLHS(mt * 2 + 1, parallelismWeight*w);
			eqs.addEqLHS(ms * 2 + 1, -parallelismWeight*w);
			eqs.addEqLHS(2 * n + eIdx, -parallelismWeight*g[e].avgRoot.y()*w);
			eqs.addEqRHSAndFinishEquation(0);

			eqs.addIneqLHS(2 * n + eIdx, -1.0);
			eqs.addIneqRHSAndFinishEquation(0);
			
			++eIdx;
		}
	}

	//constraint points to slide only along the cluster directions
	for (size_t v = 0; v < n; ++v)
	{
		/*eqs.addExactEqLHS(2 * v, g[v].root.x());
		eqs.addExactEqLHS(2 * v + 1, g[v].root.y());
		eqs.addExactRHSAndFinishEquation(g[v].root.dot(g[v].location));*/

		Eigen::Vector2d perp(-g[v].root.y(), g[v].root.x());
		for (Eigen::Vector2d vec : {perp, g[v].root})
		{
			eqs.addIneqLHS(2 * v, vec.x());
			eqs.addIneqLHS(2 * v + 1, vec.y());
			eqs.addIneqRHSAndFinishEquation(vec.dot(g[v].location) + g[v].width / 2);

			eqs.addIneqLHS(2 * v, -vec.x());
			eqs.addIneqLHS(2 * v + 1, -vec.y());
			eqs.addIneqRHSAndFinishEquation(-vec.dot(g[v].location) + g[v].width / 2);
		}
	}

	const double regWeight = 30;
	for (int i = 0; i+1 < eIdx; ++i)
	{
		eqs.addEqLHS(2 * n + i, regWeight);
		eqs.addEqLHS(2 * n + (i + 1), -regWeight);
		eqs.addEqRHSAndFinishEquation(1.0*regWeight);
	}

	GurobiLSLinIneqSolver solver;
	auto x = solver.solveLeastSquares(eqs.eq, eqs.eqB, eqs.ineq, eqs.ineqB, eqs.exactEq, eqs.exactB, 2 * n + eIdx, 0, 10000); //10000 = max image dimension

	for (size_t i = 0; i < boost::num_vertices(g); ++i)
		g[i].location = Eigen::Vector2d(x[2*i],x[2*i+1]);

	std::cout << "Done." << std::endl;
}
