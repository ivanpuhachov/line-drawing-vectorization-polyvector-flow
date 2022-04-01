#ifndef _LSEQUATIONS_AND_CONSTRAINTS_H_
#define _LSEQUATIONS_AND_CONSTRAINTS_H_
#include "GurobiLSLinIneqSolver.h"

struct Equation
{
	typedef Eigen::Triplet<double> Triplet;
	std::vector<Triplet> lhs;
	std::vector<double> rhs;

	void add(Equation other)
	{
		auto offsetRow = [](std::vector<Triplet>& tripletsInGlobalIndices, int offset)
		{
			for (auto& tr : tripletsInGlobalIndices)
				tr = Triplet(tr.row() + offset, tr.col(), tr.value());
		};

		offsetRow(other.lhs, (int)rhs.size());
		lhs.insert(lhs.end(), other.lhs.begin(), other.lhs.end());
		rhs.insert(rhs.end(), other.rhs.begin(), other.rhs.end());
	}

	void addLHS(int index, double value)
	{
		lhs.push_back(Triplet((int)rhs.size(), index, value));
	}

	void addRHSAndFinishEquation(double value)
	{
		rhs.push_back(value);
	}
};


/**
* Least Squares Program LSProg
* Consists of soft equalities (LSeq), inequalities (Ineq), and equalities (Eq)
*/

struct LSEquationsAndConstraints
{
	enum Type {Hard, LeastSquares, Inequality};
	Equation lseq, ineq, eq;
	void clear()
	{
		lseq.lhs.clear();
		lseq.rhs.clear();
		ineq.lhs.clear();
		ineq.rhs.clear();
		eq.lhs.clear();
		eq.rhs.clear();
	}

	void add(LSEquationsAndConstraints other)
	{
		lseq.add(other.lseq);
		ineq.add(other.ineq);
		eq.add(other.eq);
	}
	
	void addLHS(int index, double value, Type type)
	{
		switch (type)
		{
		case Hard:
			addEqLHS(index, value);
			break;
		case LeastSquares:
			addLSeqLHS(index, value);
			break;
		case Inequality:
			addIneqLHS(index, value);
			break;
		}
	}

	void addRHSAndFinishEquation(double value, Type type)
	{
		switch (type)
		{
		case Hard:
			addEqRHSAndFinishEquation(value);
			break;
		case LeastSquares:
			addLSeqRHSAndFinishEquation(value);
			break;
		case Inequality:
			addIneqRHSAndFinishEquation(value);
			break;
		}
	}
	std::vector<double> solve(int nVars, double variableScale);
private:
	void addLSeqLHS(int index, double value) { lseq.addLHS(index, value); }
	void addIneqLHS(int index, double value) { ineq.addLHS(index, value); }
	void addEqLHS(int index, double value) { eq.addLHS(index, value); }

	void addLSeqRHSAndFinishEquation(double value) { lseq.addRHSAndFinishEquation(value); }
	void addIneqRHSAndFinishEquation(double value) { ineq.addRHSAndFinishEquation(value); }
	void addEqRHSAndFinishEquation(double value) { eq.addRHSAndFinishEquation(value); }
	GurobiLSLinIneqSolver solver;
};

#endif