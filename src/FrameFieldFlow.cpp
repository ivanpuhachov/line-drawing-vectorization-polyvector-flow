#include "stdafx.h"
#include "FrameFieldFlow.h"
#include <complex>
#include <Eigen/Sparse>
#include <math.h>
#include "Simplify.h"
#include "Resample.h"
#include "FrameFieldFlowAux.h"
#include <iostream>
#include <fstream>

double FrameFieldFlow::sum_gaussian_at_pt(Eigen::Vector2d pt, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const
{
    const int N = initialPoly.size();
    double A = 1.0;
    auto y = pt[1];
    auto x = pt[0];
    double sum = 0;
    for(int i = 0; i < N; i++) {
        auto sigma = reebGraph[chain[i]].solvedWidth / 6;
        auto x0 = initialPoly[i].x();
        auto y0 = initialPoly[i].y();
        if (Eigen::Vector2d(x - x0, y - y0).squaredNorm() > 16 * sigma * sigma)
            continue;

        auto one_gaussian = A * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
        sum += one_gaussian;

        // Check for long segments of the polyline. If longer than sigma, sample more gaussians along the segment.
        if (i < N - 1) {
            auto total_dist = (initialPoly[i+1] - initialPoly[i]).norm();
            if (total_dist > sigma) {
                auto unit_dir = (initialPoly[i+1] - initialPoly[i]) / total_dist;
                auto steps = ceil(total_dist/sigma);
                auto step_dist = total_dist / steps;
                for(int j = 1; j < steps; j++) {
                    auto x0 = (initialPoly[i] + (j * step_dist) * unit_dir).x();
                    auto y0 = (initialPoly[i] + (j * step_dist) * unit_dir).y();
                    auto one_gaussian = A * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
                    sum += one_gaussian;
                }
            }
        }
    }
    return sum;
}

Eigen::Vector2d FrameFieldFlow::d_sum_gaussian_at_pt_wr2pt(Eigen::Vector2d pt, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const
{
    Eigen::Vector2d derivative;
    derivative << 0, 0;
    double gaussian_sum = FrameFieldFlow::sum_gaussian_at_pt(pt, initialPoly, chain);
    if (gaussian_sum >= 1.0 || gaussian_sum <= 0) {
        return derivative;
    }
    const int N = initialPoly.size();
    double A = 1.0;
    auto x = pt[0];
    auto y = pt[1];
    for(int i = 0; i < N; i++) {
        auto sigma = reebGraph[chain[i]].solvedWidth / 6;
        auto x0 = initialPoly[i].x();
        auto y0 = initialPoly[i].y();

        if (Eigen::Vector2d(x - x0, y - y0).squaredNorm() > 16 * sigma * sigma)
            continue;

        auto d_gaussian_dx = - A / (2 * pow(sigma, 2)) * (2 * (x - x0)) * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
        auto d_gaussian_dy = - A / (2 * pow(sigma, 2)) * (2 * (y - y0)) * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
        derivative(0) += d_gaussian_dx;
        derivative(1) += d_gaussian_dy;

        // Check for long segments of the polyline. If longer than sigma, sample more gaussians along the segment.
        if (i < N - 1) {
            auto total_dist = (initialPoly[i+1] - initialPoly[i]).norm();
            if (total_dist > sigma) {
                auto unit_dir = (initialPoly[i+1] - initialPoly[i]) / total_dist;
                auto steps = ceil(total_dist/sigma);
                auto step_dist = total_dist / steps;
                for(int j = 1; j < steps; j++) {
                    auto x0 = (initialPoly[i] + (j * step_dist) * unit_dir).x();
                    auto y0 = (initialPoly[i] + (j * step_dist) * unit_dir).y();
                    auto d_gaussian_dx = - A / (2 * pow(sigma, 2)) * (2 * (x - x0)) * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
                    auto d_gaussian_dy = - A / (2 * pow(sigma, 2)) * (2 * (y - y0)) * exp( - ((pow((x - x0), 2) + pow((y - y0), 2)) / (2 * pow(sigma, 2))));
                    derivative(0) += d_gaussian_dx;
                    derivative(1) += d_gaussian_dy;
                }
            }
        }
    }
    return derivative;
}

std::vector<MyPolyline> FrameFieldFlow::flow(const MyPolyline& initialPoly, std::vector<double> dist, std::vector<vertex_descriptor> chain)
{
    int NIter = 1500;
    return FrameFieldFlow::flow(initialPoly, dist, NIter, chain);
}


std::vector<MyPolyline> FrameFieldFlow::flow(const MyPolyline& initialPoly, std::vector<double> dist, const int NIter, std::vector<vertex_descriptor> chain)
{
	using namespace std;
    using namespace std::complex_literals;
    //std::cout<<"chain len "<< chain.size() <<std::endl;
    //std::cout<<"mypoly len "<<initialPoly.size()<<std::endl;
    //std::cout<<chain[0]<<std::endl;
    //std::cout<<initialPoly[0]<<std::endl;
    //std::cout<<reebGraph[chain[0]].location<<std::endl;

    bool dynamic_alpha = false;
	//double wReg = 0.1;

	const int old_N = initialPoly.size();
    //std::cout << "nb of vertices : " << old_N << std::endl;
	if (old_N <= 2)
		return { initialPoly };

    //std::cout << "Initial polyline analysis ." << std::endl;
    //double len = total_edge_length(initialPoly);
    int len = dist.size();
    MyPolyline resamp0 = resample(initialPoly, len);
    //std::cout << "Resampled polyline analysis ." << std::endl;
    //analysis(resamp0);

	//std::vector<MyPolyline> result = { initialPoly };
    //const int N = old_N;

    std::vector<MyPolyline> result = { resamp0 };
    const int N = resamp0.size();
    double avgEdgeLength = total_edge_length(resamp0) / (N-1);

    bool gi_max_went_down = false;
    int cutoff_iter_since_resample = 3;
    int iter_since_resample = 0;

    //std::cout << "nb of vertices : " << N << std::endl;

    //if (N==4) NIter = 3;//176;
    //if (N==108) NIter  = 181;//40;//340; //181;
    //if (N==310) NIter = 179;

    double min_g = 100000000.0;
    double min_q = 100000000.0;
    double min_e = 100000000.0;
    double previous_e = 544.201;
    double previous_gi_max;
    int min_g_idx;
    int min_q_idx;
    int min_e_idx;

    Eigen::VectorXd previous_g;
    int i;
    for (i = 0; i < NIter; ++i)
    {
	    double tmp_e = fullEnergy(result.back());
        double tmp_g = fullEnergy_g(result.back());
        double tmp_q = fullEnergy_q(result.back());

        if (tmp_e < min_e) {
            min_e = tmp_e;
            min_e_idx = i;}
        if (tmp_g < min_g) {
            min_g = tmp_g;
            min_g_idx = i;}
        if (tmp_q < min_q) {
            min_q = tmp_q;
            min_q_idx = i;}



	    //std::cout << i << std::endl;
        //std::cout << " ENERGY         : "<< tmp_e <<std::endl;
        //std::cout << " ENERGY (only g): "<< tmp_g <<std::endl;
        //std::cout << " ENERGY (only q): "<< tmp_q <<std::endl;
        //std::cout << " ENERGY DELTA   : "<< tmp_e - previous_e <<std::endl;
        //outfile<<tmp_e<< " " << tmp_g << " " << tmp_q << " " << tmp_e - previous_e << std::endl;


        previous_e = tmp_e;

	    if (i==NIter-1)
        {
	       //std::cout<<"LAST ITERATION , MAX INDEX : "<<max_edge_index(result.back())<<std::endl;
        }
        //std::cout << "N : " << N << std::endl;
        auto& curPoly = result.back();
        //analysis(curPoly);
        if (max_edge_length(result.back()) > 2.0 * avgEdgeLength) {
            //std::cout << "RESAMPLING IN ITER : " << i << " , MAX INDEX BEFORE : "<< max_edge_index(result.back())<<std::endl;
            std::cout << "RESAMPLING IN ITER : " << i << std::endl;
            //analysis(curPoly);
            // Maybe try that after
            // len *= 2;
            MyPolyline resamp_too_long = resample(result.back(), len);
            //analysis(resamp_too_long);
            curPoly = resamp_too_long;
            iter_since_resample = 0;
        }

		std::vector<Eigen::Vector2d> midpoints(N);
		std::vector<Eigen::Vector2d> tau(N);

		midpoints[0] = curPoly[0];
		for (int j = 1; j < N; ++j)
		{
			midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
			tau[j] = curPoly[j] - curPoly[j - 1];
		}
		tau[0] = tau[1];

		Eigen::VectorXcd X(N);
		for (int j = 0; j < N; ++j)
			X(j) = curPoly[j].x() - std::complex<double>(0.0, 1.0) * curPoly[j].y();

		std::vector<double> q(N);
		for (int j = 0; j < N; ++j)
			q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

		Eigen::VectorXcd tau_cmplx(N);
		for (int j = 0; j < N; ++j)
			tau_cmplx(j) = tau[j].x() - std::complex<double>(0.0, 1.0) * tau[j].y();

		std::vector<Eigen::Vector2d> midpt_gradients(N);
		for (int j = 0; j<N; j++) {
		    midpt_gradients[j] = FrameFieldFlow::d_sum_gaussian_at_pt_wr2pt(midpoints[j],initialPoly, chain);
		}

		std::vector<Eigen::Vector2d> normals(N);
		for (int j = 0; j < N; ++j)
			normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

		std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
		auto g = f(normals, midpoints, dg_dn, dg_dx);

		 //gaussians over curve correction
		for (int j = 0; j < N; j++) {
		    dg_dx[j].real() -= midpt_gradients[j];
		}

		dg_dx.back() = Eigen::Vector2d(0, 0); //fix the endpoint

		Eigen::VectorXcd gamma_prime(N);
		for (int j = 0; j < N; ++j)
			gamma_prime(j) = dg_dn[j].x() - std::complex<double>(0.0, 1.0) * dg_dn[j].y();

		//std::cout << g << std::endl;
		//std::cout<<"AAAAA"<<std::endl;
		//std::vector<double> factors = gradient_ff(curPoly);
		//if (i==0){
		//    for (double j : factors) std::cout<<j<<std::endl;
		//}
		//std::vector<double> tmp_factors;
        //if (i >= 1)
        //{
        //    factors = FrameFieldFlow::gradient_ff(curPoly);
        //    //auto min_element = *std::min_element(factors.begin(), factors.end());
        //    //auto max_element = *std::max_element(factors.begin(), factors.end());
        //    //std::cout<<"MAX GRAD (*100): " <<max_element*100<<std::endl;
        //    int window = 3;
        //    std::cout<<"N = "<<N<<std::endl;
        //    for (int j = 0; j<N; ++j)
        //    {
        //        int min_idx = j - window >=0 ? j -window: 0;
        //        int max_idx = j + window < N ? j + window: N;
        //        //std::cout<<"MIN MAX : "<<min_idx<<" "<<max_idx<<std::endl;
        //        int sum = 0;
        //        for (int k = min_idx; k < max_idx; k ++)
        //        {
        //            sum += factors[k];
        //        }
        //        sum /= (max_idx - min_idx + 1);
        //        tmp_factors.push_back(sum);
        //        //std::cout<<"    k"<<std::endl;
        //    }
        //    for (int j = 0; j<N; ++j)
        //    {
        //        tmp_factors[j] = (tmp_factors[j] * 2 + 0.3 > 1.0)? tmp_factors[j] = 1: tmp_factors[j] = tmp_factors[j] * 2 + 0.3;
        //    }
            //for(int j=0; j<N;)
        //}
		Eigen::VectorXcd G(N), F(N);
		for (int j = 0; j < N; ++j)
		{
		    //if (i>= 20)
		    //    G(j) = g(j) + std::complex<double>(0.0, 1.0) * (std::conj(tau_cmplx(j)) * gamma_prime(j) / q[j]).real() + wReg * dist[j];
            //else
            G(j) = g(j) + std::complex<double>(0.0, 1.0) * (std::conj(tau_cmplx(j)) * gamma_prime(j) / q[j]).real() + wReg;
            F(j) = dg_dx[j].dot(normals[j]);
        }

		typedef Eigen::SparseMatrix<std::complex<double>> Mat;
		typedef Eigen::Triplet<std::complex<double>> Triplet;
		Mat A(N-2,N-2);

		auto filterDistance = [](double d){
		    double maxvalue = 2;
		    double minvalue = 0.1;
		    double thr = 15;
            d < thr ? d = maxvalue : d = maxvalue*thr / d;
            d < minvalue ? d = minvalue : d = d;
            return d;
		};

		std::vector<Triplet> triplets;
		for (int j = 1; j < N-1; ++j)
		{
            //auto diagValue = (1.0 / (2 * alpha)) * (q[j] + q[j+1]) + G[j] / q[j] + G[j+1]/q[j+1] - 0.5 * std::complex<double>(0.0, 1.0) * F[j] + 0.5 * std::complex<double>(0.0, 1.0) * F[j+1];
            auto diagValue = (1.0 / (2 * alpha * filterDistance(dist[j]))) * (q[j] + q[j+1]) + G[j] / q[j] + G[j+1]/q[j+1] - 0.5 * std::complex<double>(0.0, 1.0) * F[j] + 0.5 * std::complex<double>(0.0, 1.0) * F[j+1];
            triplets.push_back(Triplet(j - 1, j - 1, diagValue));
            if (j > 1)
			{
				auto leftDiag = -G[j] / q[j] + 0.5 * std::complex<double>(0.0, 1.0) * F[j];
				triplets.push_back(Triplet(j - 1, j - 2, leftDiag));
			}

			if (j < N - 2)
			{
				auto rightDiag = -G[j + 1] / q[j + 1] - 0.5 * std::complex<double>(0.0, 1.0) * F[j + 1];
				triplets.push_back(Triplet(j - 1, j, rightDiag));
			}
		}

		Eigen::VectorXcd rhs(N - 2);
		for (int j = 1; j < N - 1; ++j)
			//rhs(j - 1) = (1.0 / (2 * alpha)) *X[j]*(q[j] + q[j + 1]);
            rhs(j - 1) = (1.0 / (2 * alpha * filterDistance(dist[j]))) *X[j]*(q[j] + q[j + 1]);

		rhs(0) -= X(0) * (-G[1] / q[1] + 0.5 * std::complex<double>(0.0, 1.0) * F[1]);
		rhs(N - 3) -= X(N - 1) * (-G[N-1] / q[N-1] - 0.5 * std::complex<double>(0.0, 1.0) * F[N-1]);

		A.setFromTriplets(triplets.begin(),triplets.end());

		//std::ofstream fOut("Amatrix.txt");
		//fOut << A;
		//fOut.close();

		//std::ofstream fOutR("rhs.txt");
		//fOutR << rhs;
		//fOutR.close();

		Eigen::SparseLU<Mat> solver(A);
		Eigen::VectorXcd Xnew = solver.solve(rhs);
		//std::cout << "Solution: " << Xnew << std::endl;

        //MyPolyline poly(initialPoly.size());
        MyPolyline poly(N);
        poly[0] = initialPoly[0];
        poly.back() = initialPoly.back();

		for (int j = 1; j < N - 1; ++j)
			poly[j] = Eigen::Vector2d(Xnew(j - 1).real(), -Xnew(j-1).imag());

		result.push_back(poly);
		int r = result.size();

		for (int j = 0; j < r; ++j){
        //    std::cout << "RESULTS LENGTH : " << result[j].size() << std::endl;
        }

		//double angle_cutoff = 15.0 + i / 75.0;
        //double max_a = max_angle(result.back());
		//if (max_a < angle_cutoff)
        //{
		//    std::cout<<"Reached minimum angle : "<<angle_cutoff<< " degrees, at iteration "<<i<<std::endl;
		//    break;
        //}
        /*
        if (false && (i >= 1))
        {
            std::cout<<"g - previous g : "<< (g - previous_g).norm() <<std::endl;
            if ((g-previous_g).norm() < 0.001)
            {
                break;
            }
        }
        previous_g = g;


        //std::cout<<fullEnergy_avg(result.back())<<" <= "<<result.back().size()<<" * "<<fullEnergy_max_edge(result.back())<<std::endl;
        double gi_max = fullEnergy_g_max_edge(result.back());

        outfile << fullEnergy_max_edge(result.back()) << " "<<gi_max<<std::endl;

        if(i >= 200 && !gi_max_went_down && previous_gi_max > gi_max){
            gi_max_went_down = true;
        }

        if(i >= 200 && previous_gi_max < gi_max && gi_max_went_down && iter_since_resample >= cutoff_iter_since_resample )
        {
            break;
        }
        //if ( i >= 1 && gi_max < 2.7){
        //    break;
        //}
        previous_gi_max = gi_max;
        iter_since_resample += 1;
        //std::cout<<fullEnergy_max_edge(result.back())<<" et "<<fullEnergy_g_max_edge(result.back())<<std::endl;
        //if  (fullEnergy_g_max_edge(result.back()) < 1.65)
        //{
        //    break;
        //}
        */
	}
	//std::cout<<"NB ITER :"<<i<<std::endl;
    //std::cout << " MINIMUM ENERGY: "<< min_e <<" at idx : "<<min_e_idx<<std::endl;
    //std::cout << " MINIMUM ENERGY avg per edge : "<< min_g / result.back().size() <<" at idx : "<<min_g_idx<<std::endl;
    //std::cout << " ENERGY avg per edge at the end: "<< fullEnergy_avg(result.back())<<" at idx : "<<min_g_idx<<std::endl;
    //std::cout << " MINIMUM ENERGY (g): "<< min_g <<" at idx : "<<min_g_idx<<std::endl;
    //std::cout << " MINIMUM ENERGY (q): "<< min_q <<" at idx : "<<min_q_idx<<std::endl;
    //std::cout << " MIN ENERGY MAX EDGE : "<< previous_gi_max <<std::endl<<"------ ------ ------"<<std::endl;

    //outfile<<std::endl;
    //outfile.close();
	return  result;//*/{ result.back() };
}

FrameFieldFlow::FrameFieldFlow(const Eigen::MatrixXcd& F1, const Eigen::MatrixXcd& F2, const G reebGraph)
	:F1_(F1),F2_(F2),alpha(1e-2),lengthRegularizer(0.2),wReg(3.0),reebGraph(reebGraph)
{
}

Eigen::VectorXd FrameFieldFlow::f(const std::vector<Eigen::Vector2d>& normals, const std::vector<Eigen::Vector2d>& pts, std::vector<Eigen::Vector2cd>& df_dnu, std::vector<Eigen::Vector2cd>& df_dx) const
{
	using namespace std;
	int N = pts.size();

	Eigen::VectorXd fVec(N);
	df_dnu.resize(N);
	df_dx.resize(N);

	for (int i = 0; i < N; ++i)
	{
		auto nu_cmplx = normals[i].x() + std::complex<double>(0.0, 1.0) * normals[i].y();
		auto tau = std::complex<double>(0.0, 1.0) * nu_cmplx;
		const auto& p = pts[i];
		//std::array<int, 2> pi = { static_cast<int>(std::round(p.y())),static_cast<int>(std::round(p.x())) };
		std::array<int, 2> pi = { static_cast<int>(std::round(p.x())),static_cast<int>(std::round(p.y())) };
		auto f1 = F1_(pi[0], pi[1]);
		auto f2 = F2_(pi[0], pi[1]);

		double f1R = f1.real(), f1I = -f1.imag(), f2R = f2.real(), f2I = -f2.imag();
		double nu1 = normals[i].x(), nu2 = normals[i].y();
		fVec[i] = std::pow(std::abs(f1 + f2*std::pow(tau,2) + std::pow(tau,4)),2);

		df_dnu[i].x() = 4 * std::pow(f2I, 2) * std::pow(nu1, 3) + 4 * std::pow(f2I, 2) * nu1 * std::pow(nu2, 2) + 20 * f2I * std::pow(nu1, 4) * nu2 +
			24 * f2I * std::pow(nu1, 2) * std::pow(nu2, 3) + 4 * f1I * f2I * nu1 + 4 * f2I * std::pow(nu2, 5) - 4 * f1R * f2I * nu2 +
			4 * std::pow(f2R, 2) * std::pow(nu1, 3) + 4 * std::pow(f2R, 2) * nu1 * std::pow(nu2, 2) + 12* f2R * std::pow(nu1, 5) + 8 * f2R * std::pow(nu1, 3) * std::pow(nu2, 2)
			- 4 * f2R * nu1 * std::pow(nu2, 4) + 4 * f1R * f2R * nu1 + 4 * f1I * f2R * nu2 + 8 * std::pow(nu1, 7) + 24 * std::pow(nu1, 5) * std::pow(nu2, 2) + 24 * std::pow(nu1, 3) * std::pow(nu2, 4)
			+ 8 * f1R * pow(nu1, 3) + 24 * f1I * pow(nu1, 2) * nu2 + 8* nu1 * pow(nu2, 6) - 24 * f1R * nu1 * pow(nu2, 2) - 8 * f1I * pow(nu2, 3)
			+ nu1 * (pow(4 * pow(nu1, 3) * nu2 + f2I * pow(nu1, 2) - 4 * nu1 * pow(nu2, 3) + 2 * f2R * nu1 * nu2 - f2I * pow(nu2, 2) + f1I, 2) +
				pow(pow(nu1, 4) - 6 * pow(nu1, 2) * pow(nu2, 2) + f2R * pow(nu1, 2) - 2 * f2I * nu1 * nu2 + pow(nu2, 4) - f2R * pow(nu2, 2) + f1R, 2));
		df_dnu[i].y() = (pow(4 * pow(nu1,3) * nu2 + f2I*pow(nu1,2) - 4 * nu1*pow(nu2,3) + 2 * f2R*nu1*nu2 - f2I*pow(nu2,2) + f1I,2) +
			pow(pow(nu1,4) - 6 * pow(nu1,2) * pow(nu2,2) + f2R*pow(nu1,2) - 2* f2I*nu1*nu2 + pow(nu2,4) - f2R*pow(nu2,2) + f1R, 2)) +
			(4 * pow(f2I, 2) * pow(nu1,2) * nu2 + 4 * pow(f2I, 2) * pow(nu2,3) + 4 * f2I*pow(nu1,5) + 24 * f2I*pow(nu1,3) * pow(nu2,2) +
				20 * f2I*nu1*pow(nu2,4) - 4 * f1R*f2I*nu1 - 4 * f1I*f2I*nu2 + 4 * pow(f2R, 2) * pow(nu1,2) * nu2 + 4 * pow(f2R, 2) * pow(nu2,3)
				+ 4 * f2R*pow(nu1,4) * nu2 - 8 * f2R*pow(nu1,2) * pow(nu2,3) + 4* f1I*f2R*nu1 - 12* f2R*pow(nu2,5) - 4* f1R*f2R*nu2
				+ 8* pow(nu1, 6) * nu2 + 24 * pow(nu1,4) * pow(nu2,3) + 8* f1I*pow(nu1,3) + 24* pow(nu1,2) * pow(nu2,5) - 24* f1R*pow(nu1,2) * nu2
				- 24* f1I*nu1*pow(nu2,2) + 8* pow(nu2,7) + 8* f1R*pow(nu2,3));

		complex<double> expr1 = sqrt(std::complex<double>(pow(nu1, 2) + 2* nu2));
		array<complex<double>, 4> df_dfield = { expr1*(2* pow(nu1,4) - 12* pow(nu1,2)* pow(nu2,2) + 2* f2R*pow(nu1,2) - 4* f2I*nu1*nu2 + 2* pow(nu2,4) - 2* f2R*pow(nu2,2) + 2* f1R),
			expr1*(8* pow(nu1,3)* nu2 + 2* f2I*pow(nu1,2) - 8* nu1*pow(nu2,3) + 4* f2R*nu1*nu2 - 2* f2I*pow(nu2,2) + 2* f1I),
			expr1*(2* pow(nu1,6) + 2* pow(nu1,4)* pow(nu2,2) + 2* f2R*pow(nu1,4) - 2* pow(nu1,2)* pow(nu2,4) + 4* f2R*pow(nu1,2)* pow(nu2,2) + 2* f1R*pow(nu1,2) + 4* f1I*nu1*nu2 - 2* pow(nu2, 6) + 2* f2R*pow(nu2,4) - 2* f1R*pow(nu2,2)),
			expr1*(4* pow(nu1,5)* nu2 + 2* f2I*pow(nu1,4) + 8* pow(nu1,3)* pow(nu2,3) + 4* f2I*pow(nu1,2)* pow(nu2,2) + 2* f1I*pow(nu1,2) + 4* nu1*pow(nu2,5) - 4* f1R*nu1*nu2+ 2* f2I*pow(nu2,4) - 2* f1I*pow(nu2,2)) };



		//dir = 0 --> df/dx, dir = 1 --> df/dy
		auto dM_dx = [](const Eigen::MatrixXcd& M, size_t i, size_t j, size_t dir)
		{
			size_t nMax = dir == 0 ? M.rows() : M.cols();
			size_t idx = dir == 0 ? i : j;
			if (idx != 0)
			{
				if (idx + 1 < nMax)
				{
					if (dir == 0)
						return 0.5 * (M(i + 1, j) - M(i - 1, j));
					else
						return 0.5 * (M(i, j + 1) - M(i, j - 1));
				}
				else //idx + 1 == nMax
				{
					if (dir == 0)
						return M(i, j) - M(i - 1, j);
					else
						return M(i, j) - M(i, j - 1);
				}
			}
			else //idx == 0
			{
				if (dir == 0)
					return M(i + 1, j) - M(i, j);
				else
					return M(i, j + 1) - M(i, j);
			}
		};

		array<complex<double>, 2> dfield1_dx = { dM_dx(F1_,pi[0],pi[1],0), dM_dx(F1_,pi[0],pi[1],1) };
		
		array<complex<double>, 2> dfield2_dx = { dM_dx(F2_,pi[0],pi[1],0), dM_dx(F2_,pi[0],pi[1],1) };
		df_dx[i] = { df_dfield[0] * dfield1_dx[0].real() - df_dfield[1] * dfield1_dx[0].imag() + df_dfield[2] * dfield2_dx[0].real() - df_dfield[3] * dfield2_dx[0].imag(),
					 df_dfield[0] * dfield1_dx[1].real() - df_dfield[1] * dfield1_dx[1].imag() + df_dfield[2] * dfield2_dx[1].real() - df_dfield[3] * dfield2_dx[1].imag() };

	}

	return fVec;
}


std::vector<double> FrameFieldFlow::gradient_ff(const MyPolyline& curPoly)
{
    using namespace std;
    int N = curPoly.size();
    std::vector<double> gradz;
    for (int i = 0; i < N; ++i)
    {
        //auto p = 0.5 * (curPoly[i] + curPoly[i + 1]);
        auto p = curPoly[i];
        std::array<int, 2> pi = { static_cast<int>(std::floor(p.x())),static_cast<int>(std::floor(p.y())) };

        //dir = 0 --> df/dx, dir = 1 --> df/dy
        auto dM_dx = [](const Eigen::MatrixXcd& M, size_t i, size_t j, size_t dir)
        {
            size_t nMax = dir == 0 ? M.rows() : M.cols();
            size_t idx = dir == 0 ? i : j;
            if (idx != 0)
            {
                if (idx + 1 < nMax)
                {
                    if (dir == 0)
                        return 0.5 * (M(i + 1, j) - M(i - 1, j));
                    else
                        return 0.5 * (M(i, j + 1) - M(i, j - 1));
                }
                else //idx + 1 == nMax
                {
                    if (dir == 0)
                        return M(i, j) - M(i - 1, j);
                    else
                        return M(i, j) - M(i, j - 1);
                }
            }
            else //idx == 0
            {
                if (dir == 0)
                    return M(i + 1, j) - M(i, j);
                else
                    return M(i, j + 1) - M(i, j);
            }
        };

        array<complex<double>, 2> dfield1_dx = { dM_dx(F1_,pi[0],pi[1],0), dM_dx(F1_,pi[0],pi[1],1) };

        array<complex<double>, 2> dfield2_dx = { dM_dx(F2_,pi[0],pi[1],0), dM_dx(F2_,pi[0],pi[1],1) };
        double grad_val = (pow(dfield1_dx[0].real(),2)) + (pow(dfield1_dx[0].imag(),2))
                        + (pow(dfield1_dx[1].real(),2)) + (pow(dfield1_dx[1].imag(),2))
                        + (pow(dfield2_dx[0].real(),2)) + (pow(dfield2_dx[0].imag(),2))
                        + (pow(dfield2_dx[1].real(),2)) + (pow(dfield2_dx[1].imag(),2));
        gradz.push_back(grad_val);
    }
    std::vector<double> factors;

    if(true) {
        int window = 1;
        //std::cout<<"N = "<<N<<std::endl;
        for (int j = 0; j < N; ++j) {
            int min_idx = j - window >= 0 ? j - window : 0;
            int max_idx = j + window < N ? j + window : N-1;
            //std::cout<<"MIN MAX : "<<min_idx<<" "<<max_idx<<std::endl;
            double sum = 0;
            for (int k = min_idx; k < max_idx; k++) {
                if (k==j) continue;
                sum += gradz[k];
            }
            sum = sum / (window * 2);
            sum = 0.65 * gradz[j] + 0.35 * (sum);

            factors.push_back(sum);
        }
        for (int j = 0; j < N; ++j) {
            factors[j] = (1 - factors[j] * 2 > 0.1) ? factors[j] = (1 - factors[j] * 2) : factors[j] = 0.1;
            //factors[j] = 1 - factors[j];
        }
    }
    return factors;
}

std::vector<MyPolyline> FrameFieldFlow::get_dg_dx(const MyPolyline& curPoly)
{
    auto N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    Eigen::VectorXcd X(N);
    for (int j = 0; j < N; ++j)
        X(j) = curPoly[j].x() - std::complex<double>(0.0, 1.0) * curPoly[j].y();

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    Eigen::VectorXcd tau_cmplx(N);
    for (int j = 0; j < N; ++j)
        tau_cmplx(j) = tau[j].x() - std::complex<double>(0.0, 1.0) * tau[j].y();


    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    auto g = f(normals, midpoints, dg_dn, dg_dx);

    // Making polylines
    std::vector<MyPolyline> gradients;
    for (int i = 0; i < N; i++)
    {
        MyPolyline tmp_poly(2);
        tmp_poly[0] = curPoly[i];

        tmp_poly[1] = curPoly[i];
        //std::cout<<"XY TEST "<<tmp_poly[1].x() <<" "<< tmp_poly[1].x()<<std::endl;

        // If we take either the real part or the imaginary part of dg_dx (whichever is not null)
//        if (dg_dx[i](0).imag() == 0)
//        {
//            tmp_poly[1].x() += dg_dx[i](0).real() * 1.0;
//            tmp_poly[1].y() += dg_dx[i](1).real() * 1.0;
//        }
//        else
//        {
//            tmp_poly[1].x() += dg_dx[i](0).imag() * 1.0;
//            tmp_poly[1].y() += dg_dx[i](1).imag() * 1.0;
//        }

        // If we take always the real part of dg_dx:
        tmp_poly[1].x() += dg_dx[i](0).real() * g[i] * 2.0;
        tmp_poly[1].y() += dg_dx[i](1).real() * g[i] * 2.0;

        //std::cout<<" something "<<dg_dx[i](0).real()<<" "<<dg_dx[i](0).imag()<<" "<<dg_dx[i](1).real()<<" "<<dg_dx[i](1).imag()<<std::endl;

        gradients.push_back(tmp_poly);
    }
    return gradients;
}


double FrameFieldFlow::energyPerEdge(Eigen::Vector2d p1, Eigen::Vector2d p2) const
{
    Eigen::Vector2d midpoint;
    Eigen::Vector2d tau;

    midpoint = 0.5*(p1+p2);
    tau = p2 - p1;

    double q;
    q = sqrt(tau.squaredNorm() + lengthRegularizer);

    Eigen::Vector2d normal;
    normal = Eigen::Vector2d(-tau.y(), tau.x()) / q;

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f({ normal }, { midpoint }, dg_dn, dg_dx);
    return (g[0] + wReg) * q;
}

double FrameFieldFlow::energyPerEdge_g(Eigen::Vector2d p1, Eigen::Vector2d p2) const
{
    Eigen::Vector2d midpoint;
    Eigen::Vector2d tau;

    midpoint = 0.5*(p1+p2);
    tau = p2 - p1;

    double q;
    q = sqrt(tau.squaredNorm() + lengthRegularizer);

    Eigen::Vector2d normal;
    normal = Eigen::Vector2d(-tau.y(), tau.x()) / q;

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f({ normal }, { midpoint }, dg_dn, dg_dx);
    //return (g[0] + wReg) * q;
    return g[0];
}

double FrameFieldFlow::energyPerEdge_q(Eigen::Vector2d p1, Eigen::Vector2d p2) const
{
    Eigen::Vector2d midpoint;
    Eigen::Vector2d tau;

    midpoint = 0.5*(p1+p2);
    tau = p2 - p1;

    double q;
    q = sqrt(tau.squaredNorm() + lengthRegularizer);

    Eigen::Vector2d normal;
    normal = Eigen::Vector2d(-tau.y(), tau.x()) / q;

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f({ normal }, { midpoint }, dg_dn, dg_dx);
    //return (g[0] + wReg) * q;
    return q;
}


double FrameFieldFlow::fullEnergy(const MyPolyline& curPoly) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    double energy = 0.0;
    for (int j=0; j<N; ++j)
    {
        energy += (g[j] + wReg) * q[j]; //i think??
    }
    return energy;
}
/*
 * This function computes FrameFieldFlow::fullEnergy + gaussian to the centerline
 */
double FrameFieldFlow::fullEnergy_white(const MyPolyline& curPoly, const MyPolyline& initialPoly, const std::vector<vertex_descriptor>& chain) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    // computing gaussian
    Eigen::VectorXd gaussian_distance_to_initial(N);
    for (int j=0; j<N; ++j){
        gaussian_distance_to_initial(j) = sum_gaussian_at_pt(midpoints[j],initialPoly, chain);
    }
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    g += gaussian_distance_to_initial;
    double energy = 0.0;
    for (int j=0; j<N; ++j)
    {
        energy += (g[j] + wReg) * q[j];
    }
    return energy;
}

double FrameFieldFlow::fullEnergy_g(const MyPolyline& curPoly) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    double energy = 0.0;
    for (int j=0; j<N; ++j)
        energy += (g[j] + wReg) * q[j]; //i think??

    double energy_g = 0.0;
    for(int j =0; j<N; ++j)
        energy_g += g[j];
    return energy_g;
}

double FrameFieldFlow::fullEnergy_q(const MyPolyline& curPoly) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    double energy = 0.0;
    for (int j=0; j<N; ++j)
        energy += (g[j] + wReg) * q[j]; //i think??

    double energy_q = 0.0;
    for(int j =0; j<N; ++j)
        energy_q += q[j];
    return energy_q;
}

double FrameFieldFlow::fullEnergy_max_edge(const MyPolyline& curPoly) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    double energy = 0.0;
    for (int j=0; j<N; ++j)
        energy += (g[j] + wReg) * q[j]; //i think??
    return energy;
}

double FrameFieldFlow::fullEnergy_g_max_edge(const MyPolyline& curPoly) const
{
    const int N = curPoly.size();
    std::vector<Eigen::Vector2d> midpoints(N);
    std::vector<Eigen::Vector2d> tau(N);

    midpoints[0] = curPoly[0];
    for (int j = 1; j < N; ++j)
    {
        midpoints[j] = 0.5 * (curPoly[j] + curPoly[j - 1]);
        tau[j] = curPoly[j] - curPoly[j - 1];
    }
    tau[0] = tau[1];

    std::vector<double> q(N);
    for (int j = 0; j < N; ++j)
        q[j] = sqrt(tau[j].squaredNorm() + lengthRegularizer);

    std::vector<Eigen::Vector2d> normals(N);
    for (int j = 0; j < N; ++j)
        normals[j] = Eigen::Vector2d(-tau[j].y(), tau[j].x()) / q[j];

    std::vector<Eigen::Vector2cd> dg_dn, dg_dx;
    Eigen::VectorXd g = f(normals, midpoints, dg_dn, dg_dx);
    double max_energy = g[0];
    for (int j=1; j<N; ++j)
        if (g[j] > max_energy)
            max_energy = g[j];
    return max_energy;
}



double FrameFieldFlow::fullEnergy_avg(const MyPolyline& curPoly) const
{
    double fe = fullEnergy(curPoly);
    int n = curPoly.size();
    return fe/n;
}

double FrameFieldFlow::fullEnergy_g_avg(const MyPolyline& curPoly) const
{
    double fe = fullEnergy_g(curPoly);
    int n = curPoly.size();
    return fe/n;
}



