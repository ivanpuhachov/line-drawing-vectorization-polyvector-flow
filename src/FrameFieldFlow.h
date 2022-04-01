#pragma once
#include "typedefs.h"

class FrameFieldFlow
{
public:
	FrameFieldFlow(const Eigen::MatrixXcd& F1, const Eigen::MatrixXcd& F2, const G reebGraph);
	double fullEnergy(const MyPolyline& initialPoly) const;
	double fullEnergy_white(const MyPolyline& curPoly, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const;
    double fullEnergy_g(const MyPolyline& initialPoly) const;
    double fullEnergy_q(const MyPolyline& initialPoly) const;
    double fullEnergy_max_edge(const MyPolyline& initialPoly) const;
    double fullEnergy_g_max_edge(const MyPolyline& initialPoly) const;
	double energyPerEdge(Eigen::Vector2d p1, Eigen::Vector2d p2) const;
    double energyPerEdge_g(Eigen::Vector2d p1, Eigen::Vector2d p2) const;
    double energyPerEdge_q(Eigen::Vector2d p1, Eigen::Vector2d p2) const;
    double fullEnergy_avg(const MyPolyline& initialPoly) const;
    double fullEnergy_g_avg(const MyPolyline& initialPoly) const;
    std::vector<MyPolyline> flow(const MyPolyline& poly, std::vector<double> dist, std::vector<vertex_descriptor> chain);
    std::vector<MyPolyline> flow(const MyPolyline& poly, std::vector<double> dist, const int NIter, std::vector<vertex_descriptor> chain);
    std::vector<double> gradient_ff(const MyPolyline& curPoly);
    std::vector<MyPolyline> get_dg_dx(const MyPolyline& curPoly);
private:
	Eigen::VectorXd f(const std::vector<Eigen::Vector2d>& normals, const std::vector<Eigen::Vector2d>& pts, std::vector<Eigen::Vector2cd>& outGradN, std::vector<Eigen::Vector2cd>& outGradX) const;
	const Eigen::MatrixXcd& F1_, &F2_;
	//double analysis(const MyPolyline& poly);
	const double alpha; //time step
	const double lengthRegularizer;
	const double wReg;
	const G reebGraph;
    double sum_gaussian_at_pt(Eigen::Vector2d pt, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const;
    Eigen::Vector2d d_sum_gaussian_at_pt_wr2pt(Eigen::Vector2d pt, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const;
};