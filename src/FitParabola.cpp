#include "FitParabola.h"
#include <array>
#include <queue>
#include <set>
#include "Eigen/Dense"
#include <Eigen/Eigenvalues>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include "opencv2/core/eigen.hpp"
#include "typedefs.h"

std::vector<std::array<int, 2>> findGeoNeighborhoodOLD(cv::Mat& mask, std::array<int, 2> seedPoint, int r)
{
	if (r < 0)
		return{};

	std::map<std::array<int, 2>, int> distance;
	std::map<std::array<int, 2>, bool> processed;

	distance[seedPoint] = 0;

	auto cur = seedPoint;

	while (cur[0] != -1)
	{
		int curDist = distance[cur];

		if (curDist < r)
		{
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
				{
					if ((i == j) && (i == 1))
						continue;

					std::array<int, 2> neighbor = { cur[0] + (i - 1),cur[1] + (j - 1) };
					cv::Rect rect(cv::Point(), mask.size());
					cv::Point pt(neighbor[1], neighbor[0]);
					if (rect.contains(pt) && (mask.at<uchar>(pt)) != 0)
					{
						//update distance to this
						int dist = curDist + std::max(abs(i - 1), abs(j - 1));
						if ((distance.find(neighbor) == distance.end()) || (distance[neighbor] > dist))
							distance[neighbor] = dist;
					}
				}
		}

		processed[cur] = true;

		int bestDist = std::numeric_limits<int>::max();
		cur = { -1,-1 };
		for (auto it : distance)
		{
			if ((!processed[it.first]) && (it.second < bestDist))
			{
				cur = it.first;
				bestDist = it.second;
			}
		}
	}

	std::vector<std::array<int, 2>> result;
	result.reserve(distance.size());
	for (auto it : distance)
		result.push_back(it.first);

	return result;
}

std::vector<std::array<int, 2>> findGeoNeighborhood(cv::Mat& mask, std::array<int, 2> seedPoint, int r)
{
	if (r < 0)
		return{};

	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS, no_property,
		property< edge_weight_t, int, property< edge_weight2_t, int > > > Graph;
	typedef graph_traits < Graph >::vertex_descriptor vertex_descriptor;
	typedef std::pair < int, int >Edge;

	std::vector<Edge> myEdges;
	int m = mask.rows, n = mask.cols;

	int V = 0;
	std::map<std::array<int, 2>, int> indices;
	std::vector<std::array<int, 2>> oneDindexTo2D;
	for (int i = 0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				indices[{i, j}] = V;
				oneDindexTo2D.push_back({ i,j });
				++V;
			}
		}

	for (int i = 0; i<m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				int myIndex = indices[{i, j}];
				if ((i < m - 1) && (mask.at<uchar>(i+1, j) != 0))
				{
					int index = indices[{i + 1, j}];
					myEdges.push_back({ myIndex,index });
				}
				if ((j < n - 1)&&(mask.at<uchar>(i, j+1) != 0))
				{
					int index = indices[{i, j + 1}];
					myEdges.push_back({ myIndex,index });
				}
				if ((i < m - 1) && (j < n - 1) && (mask.at<uchar>(i+1, j+1) != 0))
				{
					int index = indices[{i + 1, j + 1}];
					myEdges.push_back({ myIndex,index });
				}
			}
		}

	Graph g(&myEdges[0], &myEdges[0] + myEdges.size(), V);

	property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
	graph_traits < Graph >::edge_iterator e, e_end;

	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
		w[*e] = 1;

	std::vector<vertex_descriptor> p(num_vertices(g));
	std::vector<int> d(num_vertices(g));
	vertex_descriptor s = vertex(indices[seedPoint], g);

	dijkstra_shortest_paths(g, s,
		predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
		distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));

	std::vector<std::array<int, 2>> result;
	for (int i = 0; i < V; ++i)
	{
		if (d[i] <= r)
			result.push_back(oneDindexTo2D[i]);
	}
	return result;
}

bool fitParabola(const cv::Mat & bwImg, float strokeRadius, const cv::Mat & mask, double x0, double y0, bool onlyReliable, Eigen::Vector2d& outCenterPt, Eigen::Vector2d& outCenterAxis, std::array<double,2> & outD, double & outResidual, int & outNotNull)
{
//#define MY_DEBUG 1

	using namespace cv;

	cv::Rect bigBounds(Point(), bwImg.size());

	outResidual = std::numeric_limits<double>::max();

	std::array<int,2> p0i = { static_cast<int>(std::round(y0)),static_cast<int>(std::round(x0)) };
	int windowSize = 2 * std::ceil(strokeRadius) + 1;
	if (windowSize % 2 == 0)
		windowSize++;

	int step = (windowSize - 1) / 2;
	int m = bwImg.rows;
	int n = bwImg.cols;

	int iMin = std::max(0, p0i[0] - step);
	int iMax = std::min(m-1, p0i[0] + step);
	int jMin = std::max(0, p0i[1] - step);
	int jMax = std::min(n-1, p0i[1] + step);

	auto myIndex = p0i;
	if (p0i[0] > step)
		myIndex[0] = step;
	if (p0i[1] > step)
		myIndex[1] = step;

	Rect rect = Rect(jMin, iMin, jMax - jMin + 1, iMax - iMin + 1);
	Mat bwCut(bwImg, rect);
	assert(bwCut.at<uchar>(myIndex[0], myIndex[1]) == bwImg.at<uchar>(p0i[0], p0i[1]));

	Mat nonNull(mask, rect);

#ifdef MY_DEBUG
	Eigen::MatrixXd notNullEigen;
	cv2eigen(nonNull, notNullEigen);
	std::cout << notNullEigen << std::endl;
#endif

	std::vector<std::array<int, 2>> geodesicNeighborhood = findGeoNeighborhood(nonNull, myIndex, step);
	int N = geodesicNeighborhood.size();
	if (N < 6)
		return false;

	Eigen::VectorXd x(N), y(N), b(N);

	for (int idx = 0; idx < N; ++idx)
	{
		auto p = geodesicNeighborhood[idx];
		y[idx] = static_cast<double>(iMin + p[0]);
		x[idx] = static_cast<double>(jMin + p[1]);
		b[idx] = static_cast<double>(bwImg.at<uchar>(iMin + p[0], jMin + p[1]));
	}

#ifdef MY_DEBUG
	std::cout << "x = " << x << std::endl;
	std::cout << "y = " << y << std::endl;
#endif

	Eigen::MatrixXd A(N, 6);
	A.col(0) = x.array().square();
	A.col(1) = y.array().square();
	A.col(2) = x.array()*y.array();
	A.col(3) = x;
	A.col(4) = y;
	A.col(5) = Eigen::VectorXd::Ones(N);
	Eigen::VectorXd coeffs = (A.transpose() * A).ldlt().solve(A.transpose() * b);

#ifdef MY_DEBUG
	std::cout << "Coeffs: " << coeffs << std::endl;
#endif

	double res = (A*coeffs - b).norm();
	Eigen::MatrixXd coefM(2,2);
	coefM << coeffs(0), coeffs(2) / 2, coeffs(2) / 2, coeffs(1);

	Eigen::EigenSolver<Eigen::Matrix2d> es(coefM);
	std::array<std::complex<double>, 2> newD = { es.eigenvalues()[0],es.eigenvalues()[1] };

	assert((fabs(newD[0].imag()) < 1e-16) && (fabs(newD[1].imag()) < 1e-16)); //should be real

	Eigen::Vector2d newCenterAxis = newD[0].real() < newD[1].real() ? es.eigenvectors().col(1).real() : es.eigenvectors().col(0).real();

	Eigen::Vector2d rhs;
	rhs << coeffs(3), coeffs(4);
	Eigen::Vector2d newCenterPt = -(coefM + coefM.transpose()).inverse()*rhs;
	Eigen::Vector2d seedPtVec(x0, y0);
	newCenterPt += (seedPtVec - newCenterPt).transpose()*newCenterAxis*newCenterAxis;

#ifdef MY_DEBUG
	bwCut.convertTo(bwCut, CV_32F,1.0/255);
	Mat bwCutUpsampled;
	const double scale = 100;
	resize(bwCut, bwCutUpsampled, Size(), scale, scale, INTER_NEAREST);
	cvtColor(bwCutUpsampled, bwCutUpsampled, CV_GRAY2RGB);
	circle(bwCutUpsampled, Point((newCenterPt.x()-jMin+0.5)*scale, (newCenterPt.y()-iMin+0.5)*scale), 10, Scalar(1, 0, 0),3);
	imshow("Display window", bwCutUpsampled);
	cv::waitKey(0);
#endif
	
	double minEig = std::min(newD[0].real(), newD[1].real()), maxEig = std::max(newD[0].real(), newD[1].real());
	bool reliable = (minEig < 0) && (fabs(maxEig/minEig)<0.1);



	Point cvNewCenterPt(std::round(newCenterPt.x()), std::round(newCenterPt.y()));

	if ((reliable || !onlyReliable) && bigBounds.contains(cvNewCenterPt) && (mask.at<uchar>(cvNewCenterPt) != 0))
	{
		outResidual = res;
		outNotNull = N;
		outCenterPt = newCenterPt;
		outCenterAxis = newCenterAxis;
		outD = { minEig, maxEig };
		return true;
	}

	return false;
}
