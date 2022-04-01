#include "LSEquationsAndConstraints.h"

std::vector<double> LSEquationsAndConstraints::solve(int nVars, double variableScale)
{
	
	auto X_gurobi_LS = solver.solveLeastSquares(lseq.lhs, lseq.rhs, ineq.lhs, ineq.rhs, eq.lhs, eq.rhs, nVars, -variableScale, variableScale);
	
	return X_gurobi_LS;
}