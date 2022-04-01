#include "WidthSolver.h"
#include "Equations.h"
#include "LSEquationsAndConstraints.h"

void findWidths(G& g)
{
	LSEquationsAndConstraints eqs;
	const double w = 1.0;

	double maxWidth = 0.0;
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		if (maxWidth < g[v].width)
			maxWidth = g[v].width;
	}

	//form the matrix: Laplacian + initial values
	for (size_t v = 0; v < boost::num_vertices(g); ++v)	{
		G::adjacency_iterator ai, a_end;
		int deg = boost::degree(v, g);
		eqs.addLHS(v, -deg, LSEquationsAndConstraints::LeastSquares);
		for (boost::tie(ai, a_end) = boost::adjacent_vertices(v, g); ai != a_end; ++ai)
			eqs.addLHS(*ai, 1.0, LSEquationsAndConstraints::LeastSquares);
		eqs.addRHSAndFinishEquation(0.0, LSEquationsAndConstraints::LeastSquares);

		eqs.addLHS(v, w, LSEquationsAndConstraints::LeastSquares);
		eqs.addRHSAndFinishEquation(w*g[v].width, LSEquationsAndConstraints::LeastSquares);
	}

	auto newWidths = eqs.solve(boost::num_vertices(g),maxWidth);
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		g[v].solvedWidth = newWidths[v];
	}
}
