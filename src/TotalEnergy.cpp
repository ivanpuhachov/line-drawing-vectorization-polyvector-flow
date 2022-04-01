#include "stdafx.h"
#include "TotalEnergy.h"
#include "polynomial_energy.h"
#include "l2_regularizer.h"
#include <iostream>
#include <fstream>
#include "Params.h"
TotalEnergy::TotalEnergy(cv::Mat & bwImg, const Eigen::MatrixXd & weight, const Eigen::MatrixXcd & tauNormalized, double beta, cv::Mat & mask, const Eigen::MatrixXi& indices, int nnz)
	:bwImg(bwImg),weight(weight),tauNormalized(tauNormalized),beta(beta),mask(mask),noisy(noisy),indices(indices),nnz(nnz)
{
	m = bwImg.rows;
	n = bwImg.cols;

	alpha = FRAME_FIELD_REGULARIZER_WEIGHT;
	smartWeights = Eigen::MatrixXd(m,n);
	smartWeights.setOnes();
	for (int i=1; i<m-1; ++i)
		for (int j = 1; j < n-1; ++j)
		{
			smartWeights(i, j) = 0.25*((1-weight(i - 1, j)) + (1-weight(i + 1, j)) + (1-weight(i, j - 1)) + (1-weight(i, j + 1)));
		}

	onesMatrix = Eigen::MatrixXd(m, n);
	onesMatrix.setOnes();

	g = tauNormalized * std::complex<double>(0, 1);
}

double TotalEnergy::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	Eigen::VectorXcd g1, g2, g3;
	double e1, e2, e3;
	Eigen::VectorXcd x_complex = x.head(x.size() / 2) + std::complex<double>(0, 1) * x.tail(x.size() / 2);
	std::vector<double> energiesOut;
	std::tie(e1, g1) = polynomial_energy(x_complex, weight, tauNormalized, mask, indices, energiesOut, true);
	//auto [test, A, b] = polynomial_energy_matrix(x_complex, weight, tauNormalized, mask, indices, energiesOut, true);

	
	//ff << "Are = [" << A.real() << "];" << std::endl;
	//ff << "Aim = [" << A.imag() << "];" << std::endl;
	//ff << "bRe = [" << b.real() << "];" << std::endl;
	//ff << "bIm = [" << b.imag() << "];" << std::endl;
	//ff << "A = Are + 1i*Aim;" << std::endl;
	//ff << "b = bRe + 1i*bIm;" << std::endl;
	//ff << "Xre = [" << x_complex.real() << "];" << std::endl;
	//ff << "Xim = [" << x_complex.imag() << "];" << std::endl;
	//ff << "X = Xre + 1i*Xim";


	//double doubleTest = test - e1;
	//Eigen::VectorXcd gradNew = (2 * A * x_complex + 2 * b.conjugate());

	//double testGrad = (g1 - gradNew).norm();
	std::tie(e2, g2) = polynomial_energy(x_complex, onesMatrix, g, mask, indices, energiesOut, true);
	int c = cv::countNonZero(mask);
	std::tie(e3, g3) = l2_regularizer(x_complex, smartWeights, mask, indices, true);


	Eigen::VectorXcd x0 = x_complex.head(x_complex.size() / 2), x2 = x_complex.tail(x_complex.size() / 2);

	Eigen::VectorXd smartWeightsVector(x_complex.size() / 2);

	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			if (mask.at<uchar>(i, j) == 0)
				continue;
			int idx = indices(i, j);
			smartWeightsVector(idx) = smartWeights(i, j);
		}
	}

	//Eigen::SparseMatrix<double> L = laplacian_matrix(mask, indices, smartWeightsVector);
	//std::complex<double> e2_matrix = (x0.adjoint() * L * x0)(0) + (x2.adjoint()* L * x2)(0);
	//Eigen::VectorXcd gradLapl (x0.size()*2);
	//gradLapl << 2.0 * L * x0, 2.0 * L * x2;
	//double testGrad2 = (g3 - gradLapl).norm();
	
	double totalEnergy = e1 + alpha*e2 + beta*e3;
	Eigen::VectorXcd grad_complex = g1 + alpha*g2 + beta*g3;
	grad.setZero();
	grad << grad_complex.real() , grad_complex.imag();

	return totalEnergy;
}

double TotalEnergy::energyOnly(const Eigen::VectorXd & x)
{
	Eigen::VectorXcd g1, g2, g3;
	double e1, e2, e3;
	Eigen::VectorXcd x_complex = x.head(x.size() / 2) + std::complex<double>(0, 1)*x.tail(x.size() / 2);
	std::vector<double> energiesOut;
	std::tie(e1, g1) = polynomial_energy(x_complex, weight, tauNormalized, mask, indices, energiesOut, false);

	std::tie(e2, g2) = polynomial_energy(x_complex, onesMatrix, g, mask, indices, energiesOut, false);
	std::tie(e3, g3) = l2_regularizer(x_complex, smartWeights, mask, indices, false);

	double totalEnergy = e1 + alpha*e2 + beta*e3;
	return totalEnergy;
}