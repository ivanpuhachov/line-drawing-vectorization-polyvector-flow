#include "stdafx.h"
#ifdef QT_VERSION
#include <QtCore/QCoreApplication>
#include "mainwindow.hpp"
#include <QElapsedTimer>
#include "mainwindow.hpp"
#endif
#ifndef _WIN32
#include <sys/time.h>
#include <unistd.h>
#endif

#include <random>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include "Eigen/Dense"
#include "opencv2/core/eigen.hpp"
#include "Eigen/Core"
#include <fstream>
#include "Optimizer.h"
#include <array>
#include "FindRoots.h"
#include "traceAuto.h"
#include "chopFakeEnds.h"
#include "simple_svg_1.0.0.hpp"
#include "Simplify.h"
#include "AlmostReebGraph.h"
#include "ContractLoops.h"
#include "RemoveShortBranches.h"
#include "SplitEmUp.h"
#include "polynomial_energy.h"
#include "findSingularities.h"
#include "Params.h"
#include "ContractDeg2.h"
#include "Smooth.h"
#include "ContractLoops2.h"
#include "FrameFieldFlow.h"
#include "WidthSolver.h"
#include "TopologyOptimizer.h"
#include "FindAllPaths.h"
#include "ExportSVG.h"
#include "ShortestPathConstrained.h"
#include "distanceGraph.h"
//#include "ChainUtils.h"
#include "SelectCurves.h"

cv::Mat bwImg, origMask;
int m, n;
Eigen::MatrixXcd g, tau, tauTimesGmag;
Eigen::MatrixXd gMag, weight;
Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

void calculateGradient()
{
	using namespace cv;
	Mat grad_x, grad_y;
	const double scale = 1.0;
	const double delta = 0;

	Sobel(bwImg, grad_x, CV_32F, 1, 0, 3, scale, delta, BORDER_DEFAULT);
	Sobel(bwImg, grad_y, CV_32F, 0, 1, 3, scale, delta, BORDER_DEFAULT);

	Mat planes[] = { grad_x, grad_y };
	Mat cvG;
	merge(planes, 2, cvG);
	cv2eigen(cvG, g);

	tauTimesGmag = g * std::complex<double>(0.0, 1.0);
	gMag = tauTimesGmag.cwiseAbs();

	double maxGradMag = gMag.maxCoeff();
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j) / maxGradMag) < 0.1)
			{
				gMag(i, j) = 0;
				tauTimesGmag(i, j) = 0;
			}
		}

	Eigen::MatrixXd gMagNoZeros = gMag;

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j)) < 1e-10)
				gMagNoZeros(i, j) = 1;
		}

	tau = tauTimesGmag.array() / gMagNoZeros.array();
}

void calculateWeight()
{
	using namespace cv;

	Eigen::MatrixXcd eigTauTimesGmag2 = tauTimesGmag.array().pow(2);

	//calculate laplacian
	//todo: get rid of so many copies back and forth
	Mat eigTauTimesGmag2Re, eigTauTimesGmag2Im;

	Eigen::MatrixXd eigTauTimesGmag2Real = eigTauTimesGmag2.real(), eigTauTimesGmag2Imag = eigTauTimesGmag2.imag();
	eigen2cv(eigTauTimesGmag2Real, eigTauTimesGmag2Re);
	eigen2cv(eigTauTimesGmag2Imag, eigTauTimesGmag2Im);

	Mat cv_gMagPow2;
	Eigen::MatrixXd gMagPow2 = gMag.array().pow(2);
	eigen2cv(gMagPow2, cv_gMagPow2);

	Mat Lx, Ly;
	Mat kernel;
	kernel = Mat::ones(3, 3, CV_64F);
	kernel.at<double>(1, 1) = 0;
	filter2D(eigTauTimesGmag2Re, Lx, -1, kernel);
	filter2D(eigTauTimesGmag2Im, Ly, -1, kernel);

	Eigen::MatrixXd Lx_eig, Ly_eig;
	cv2eigen(Lx, Lx_eig);
	cv2eigen(Ly, Ly_eig);

	Eigen::MatrixXcd mse = Lx_eig + std::complex<double>(0, 1) * Ly_eig;
	Eigen::MatrixXd mseNorm = mse.cwiseAbs();
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			if (mseNorm(i, j) < 1e-10)
				mseNorm(i, j) = 1;

	mse = mse.array() / mseNorm.array();

	Eigen::MatrixXcd tau2 = tau.array().pow(2);
	mse = mse - tau2;

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			if (fabs(gMag(i, j)) < 1e-10)
				mse(i, j) = 0;

	weight = mse.cwiseAbs();

	weight = Eigen::MatrixXd::Ones(m, n) - weight / weight.maxCoeff();

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (fabs(gMag(i, j)) < 1e-10)
			{
				weight(i, j) = 0;
			}
		}
}

Eigen::MatrixXi calculateIndices(const cv::Mat& mask, int& nnz)
{
	Eigen::MatrixXi indices = Eigen::MatrixXi(m, n);
	indices.setConstant(-1);
	int idx = 0;
	int tt = mask.type();
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < m; ++i)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				indices(i, j) = idx;
				idx++;
			}
		}
	nnz = idx;
	return indices;
}

void computeAllGeodesicDistances(const cv::Mat& mask, const Eigen::MatrixXi& indices, int nnz)
{
	std::cout << "Computing all geodesic distances... ";
	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS, no_property,
		property< edge_weight_t, int, property< edge_weight2_t, int > > > Graph;
	typedef std::pair < int, int >Edge;

	int V = nnz;
	std::vector<Edge> myEdges;
	int m = mask.rows, n = mask.cols;
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (mask.at<uchar>(i, j) != 0)
			{
				for (int sign : {-1, 1}) //leftRight
				{
					for (int dir = 0; dir < 2; dir++) //horizontal or vertical
					{
						std::pair<int, int> neigh;
						if (useNeighbor(i, j, m, n, (dir == 1), sign == -1, mask, neigh))
						{
							int idx1 = indices(i, j);
							int idx2 = indices(neigh.first, neigh.second);
							if (idx1 > idx2)
								myEdges.push_back({ idx1,idx2 });
						}
					}
				}
			}
		}

	Graph g(&myEdges[0], &myEdges[0] + myEdges.size(), V);

	property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
	graph_traits < Graph >::edge_iterator e, e_end;

	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
		w[*e] = 1;

	std::vector < int >d(V, (std::numeric_limits < int >::max)());
	std::vector<std::vector<int>> D(V, std::vector<int>(V));
	johnson_all_pairs_shortest_paths(g, D, distance_map(&d[0]));
	std::cout << "done." << std::endl;
}

void repairMask(cv::Mat& origMask)
{
	std::vector<std::pair<int, int>> newGuys;
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			if (origMask.at<uchar>(i, j) == 0)
			{
				int nn = 0;
				for (int i1 = -1; i1 < 2; ++i1)
					for (int j1 = -1; j1 < 2; ++j1)
					{
						if ((i1 == j1) || (i1 + i < 0) || (i1 + i >= m) || (j1 + j < 0) || (j1 + j >= n))
							continue;

						if (origMask.at<uchar>(i1 + i, j1 + j) != 0)
							nn++;
					}
				if (nn >= 5)
				{
					newGuys.push_back({ i, j });
				}
			}
		}

	for (auto p : newGuys)
		origMask.at<uchar>(p.first, p.second) = 255;
}

std::string type2str(int type) {
	std::string r;

	uchar depth = type & CV_MAT_DEPTH_MASK;
	uchar chans = 1 + (type >> CV_CN_SHIFT);

	switch (depth) {
	case CV_8U:  r = "8U"; break;
	case CV_8S:  r = "8S"; break;
	case CV_16U: r = "16U"; break;
	case CV_16S: r = "16S"; break;
	case CV_32S: r = "32S"; break;
	case CV_32F: r = "32F"; break;
	case CV_64F: r = "64F"; break;
	default:     r = "User"; break;
	}

	r += "C";
	r += (chans + '0');

	return r;
}

int main(int argc, char* argv[])
{
#ifdef QT_VERSION
	QApplication a(argc, argv);
#endif
	using namespace cv;
	Mat image;

    std::string filename = "../inputs/example.png";
    std::string filename_points = "../inputs/example.pts";
    int flowIterations = 1500;
    std::string svg_filename = filename + ".svg";
    bool drawMenu = true;

	if (argc < 2)
	{
		std::cout << "Usage: " << std::endl;
		std::cout << "polyvector_thing.exe FILE.png FILE.png.pts OUTFILE.svg" << std::endl;
		std::cout << "Using default \n \n" << std::endl;
//		return 0;
	}
	if (argc == 2)
	{
		filename = argv[1];
		filename_points = filename + ".pts";
        svg_filename = filename + ".svg";
	}
	else if (argc == 3)
	{
		filename = argv[1];
		filename_points = argv[2];
        svg_filename = filename + ".svg";
	}
	else if (argc == 4)
    {
        filename = argv[1];
        filename_points = argv[2];
        svg_filename = argv[3];
        drawMenu = false;
    }
	else if (argc > 4)
    {
        filename = argv[1];
        filename_points = argv[2];
        svg_filename = argv[3];
        flowIterations = std::atoi(argv[4]);
        drawMenu = false;
    }

	std::cout << "\n\n \t\t\t ---- " << filename << " ---- \n";
	std::cout << "Loading " << filename << "... ";
	image = imread(filename, IMREAD_COLOR);

	if (!image.data) // Check for invalid input
	{
		std::cout << "Could not open or find the image" << std::endl;
		return -1;
	}
	std::cout << "ok" << std::endl;

	cvtColor(image, bwImg, cv::COLOR_BGR2GRAY);
	bwImg = Scalar(255) - bwImg;
	m = bwImg.rows; n = bwImg.cols; //m: height, n: width

#ifdef _WIN32
	QElapsedTimer timer;
	timer.start();
#else
	struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL);
#endif

	//fill in the mask
	cv::Mat actualMask;
	threshold(bwImg, actualMask, BACKGROUND_FOREGROUND_THRESHOLD, 255, THRESH_BINARY);

	for (int i = 0; i < 3; ++i)
		repairMask(actualMask);

	origMask = actualMask;

	MainWindow mw;
	if (drawMenu){
        mw.setImage(QString::fromStdString(filename), actualMask);
        mw.show();
	}

	calculateGradient();
	calculateWeight();
	int nnz = 0;
	auto indices = calculateIndices(origMask, nnz);

	std::vector<double> ws;
	double maxRes = 0;

	std::vector<std::vector<Eigen::Vector2d>> centersForI(m), axiForI(m);
	std::vector<std::vector<double>> resForI(m);
	std::vector<std::vector<int>> notNullsForI(m);
	std::map<std::array<int, 2>, CenterFit> fits;

	std::cout << "Optimizing...";

	//QElapsedTimer timer1;
	//timer1.start();
//	Eigen::VectorXcd X = optimize(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
	//std::cout << "Time: " << timer1.elapsed() << std::endl;

	//QElapsedTimer timer2;
	//timer2.start();
	Eigen::VectorXcd X  = optimizeByLinearSolve(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
	//std::cout << "Time 2: " << timer2.elapsed() << std::endl;
	//std::cout << "DIFF NORM: " << (X - X2).norm() << std::endl;

	std::cout << "Finding roots.. ";
	auto roots = findRoots(X, origMask);

	auto singularities = findSingularities(roots, X, indices, origMask);
	bool improved;
	int totalNSingularities = 0;
	do
	{
		int origSingularityCount = singularities.size();

		bool somethingNew = false;
		for (auto s : singularities)
		{
			if (weight(s[0], s[1]) > 1e-5)
			{
				somethingNew = true;
				weight(s[0], s[1]) = 0;
				totalNSingularities++;
			}
		}
		if (!somethingNew)
			break;

		//X = optimize(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
		X = optimizeByLinearSolve(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, origMask, indices);
		roots = findRoots(X, origMask);
		singularities = findSingularities(roots, X, indices, origMask);

		std::cout << "done (" << origSingularityCount - (int)singularities.size() << " singularities removed" << std::endl;
		improved = origSingularityCount - singularities.size() > 0;
	} while (improved);

	mw.setRoots(roots);

	//now trace
	std::cout << "Tracing... ";
	std::map<std::array<int, 2>, std::vector<PixelInfo>> pixelInfo;
	std::vector<std::array<bool, 2>> endedWithASingularity;
	auto polys = traceAll(bwImg, actualMask, actualMask, roots, X, indices, pixelInfo, endedWithASingularity);

	mw.setPolys(polys);
	//return a.exec();
	auto singularitiesWithMask = findSingularities(roots, X, indices, origMask);
	auto reebGraph = computeAlmostReebGraph(origMask, roots, polys, pixelInfo, singularities, indices, X, endedWithASingularity);
	findWidths(reebGraph);
	
	mw.setGraph("Orig graph", reebGraph);

	contractSingularityBranches(reebGraph);
	simpleThresholds(reebGraph);
	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(reebGraph); eit != eend; ++eit)
	{
		reebGraph[*eit].weight = 1;
	}
	contractLoops(reebGraph, origMask, polys);

	mw.setGraph("Sing contracted", reebGraph);


	// find intersections/endpoints
	std::vector<Point2d> pts;
	std::vector<int> pts_types;

	std::ifstream fPts(filename_points);
	int nPts = 0;
	if (fPts.good())
	{
		fPts >> nPts;
		for (int i = 0; i < nPts; ++i)
		{
			double x, y;
			double point_type;
			keypointType kt = endpoint;
			fPts >> x >> y >> point_type;
			pts.emplace_back(x, y);
			pts_types.push_back(point_type);
		}
		fPts.close();
	}
	else
	{
		std::cout << "ERROR: pts file not found. Skipping the rest of the pipeline" << std::endl;
		a.exec();
	}



	//make sure the edge weights are correct
	for (auto [eit, eend] = boost::edges(reebGraph); eit != eend; ++eit)
	{
		reebGraph[*eit].weight = (reebGraph[eit->m_target].location - reebGraph[eit->m_source].location).norm()+1e-9;
	}

	std::map<vertex_descriptor, bool> covered;
	//addJunctions(reebGraph); //adds extra edges when strokes intersect. should be done AFTER all the width calculations

	//extend the frame field to the rest of the image
	cv::Mat fullImageMask = cv::Mat::ones(m, n, CV_8U);
	const double dilation_size = 5.0;
	cv::Mat element = getStructuringElement(cv::MORPH_ELLIPSE,
		Size(2 * dilation_size + 1, 2 * dilation_size + 1),
		Point(dilation_size, dilation_size));
	dilate(origMask, fullImageMask, element);

	int newNnz;
	auto newIndices = calculateIndices(fullImageMask, newNnz);
	auto Xextended = optimizeByLinearSolve_holdingSomeFixed(bwImg, weight, tau, FRAME_FIELD_SMOOTHNESS_WEIGHT, fullImageMask, origMask, newIndices, indices, X);


	Eigen::MatrixXcd c0(m, n), c2(m, n);

	int nonzeros = X.size() / 2;
	int idx = 0;
	for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            if (fullImageMask.at<uchar>(i, j) != 0) {
                c0(i, j) = Xextended(idx);
                c2(i, j) = Xextended(idx + newNnz);
                idx++;
            } else {
                c0(i, j) = 0;
                c2(i, j) = 0;
            }
        }
    }

    // set to true if you want to scale the framefield outside the original mask
    bool scale_ff = true;
    if(scale_ff)
    {
        Mat inverse_origMask;
        Mat dist;
        cv::bitwise_not(origMask, inverse_origMask);
        cv::distanceTransform(inverse_origMask,dist,DIST_L1,3,CV_8U);
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                if ((fullImageMask.at<uchar>(i, j) != 0) && (origMask.at<uchar>(i, j) == 0)) {
                    c0(i, j) = c0(i, j) * 2.0 * static_cast<double>(dist.at<uchar>(i, j));
                    c2(i, j) = c2(i, j) * 2.0 * static_cast<double>(dist.at<uchar>(i, j));
                    //c0(i,j) = c0(i,j) * 1.0;
                    //c2(i,j) = c2(i,j) * 1.0;
                }
            }
        }
    }

    // To view the framefield outside the original mask
    bool show_more_ff = false;
    if (show_more_ff)
    {
        roots = findRoots(Xextended, fullImageMask);
        mw.setRoots(roots);
    }

    FrameFieldFlow fff(c0, c2, reebGraph);
    std::cout << "\n\n\n FRAME FIELD " << (c0-c2).lpNorm<1>() << "\n";
    for (size_t i=0; i<5; ++i){
        std::cout << c0(i,0) << " - " << c2(i,0) << "\n";
    }
    std::vector<cv::Point2d> newlyAddedIntersections;
    std::vector<std::vector<vertex_descriptor>> finalchains, allowedChains;
    std::map<int , int> keypointToValenceMap;
    std::vector<std::pair<int, int>> allowedChainToLabeledPointPair;

    std::vector<MyPolyline> initialPolys; // rough polylines extracted from the graph
    initialPolys = optimizeTopology(reebGraph, pts, bwImg, origMask,covered,fff, roots, newlyAddedIntersections, finalchains, allowedChains, keypointToValenceMap, allowedChainToLabeledPointPair);
    auto distPerVertex = computeGraphDistancesPerPoly(reebGraph, initialPolys, allowedChains);

    mw.setCovered(covered);
    mw.setGraph("Opt graph", reebGraph);

    //    std::vector<Point2d> end_pts, sharp_pts, intersect_pts;
    std::vector<keypointType> keypoints_type;
    for (size_t i=0; i<pts.size(); ++i){
        std::cout << "Labeled point " << i << " has type " << pts_types[i] << "; computed valence is " << keypointToValenceMap[i] << "\n";
        if (pts_types[i]==1){
            if (keypointToValenceMap[i]>=1)
                keypoints_type.push_back(endpoint);
            else
                keypoints_type.push_back(error_endpoint);
        }

        if (pts_types[i]==3) {
            if (keypointToValenceMap[i]>=3)
                keypoints_type.push_back(junction);
            else
                keypoints_type.push_back(error_junction);
        }
        if (pts_types[i]==5) {
            if (keypointToValenceMap[i]>=2)
                keypoints_type.push_back(sharpcorner);
            else
                keypoints_type.push_back(error_sharpcorner);
        }
    }


    mw.setPointSet("Key points", pts, keypoints_type);

    std::vector<keypointType> newlyAdded_types(newlyAddedIntersections.size(), newpoint);
	mw.setPointSet("New intersections",newlyAddedIntersections, newlyAdded_types);

	//now run flow on them
	std::vector<MyPolyline> flowPolys(initialPolys.size()*(flowIterations+1));// = initialPolys;
	std::vector<QColor> flowColors(initialPolys.size()*(flowIterations+1)), finalColors(initialPolys.size());
	std::vector<MyPolyline> smoothedCurves(initialPolys.size()); // initialPolys after polyvector smoothing
    std::vector<MyPolyline> gradients;


    std::cout << "Starting FLOW in parallel with " << flowIterations << " iterations" << std::endl;
    #pragma omp parallel for shared(flowIterations, initialPolys, distPerVertex, flowPolys, smoothedCurves,flowColors, finalColors, c0, c2, reebGraph, finalchains) default(none)
    for (int i = 0; i < initialPolys.size(); ++i)
	{
		MyPolyline pp;
        FrameFieldFlow fff(c0, c2, reebGraph);
		for (int j = 0; j < initialPolys[i].size(); ++j)
			pp.push_back(Eigen::Vector2d(initialPolys[i][j].y(), initialPolys[i][j].x()));

        printf("CURVE : %d \n",i);
		std::vector<MyPolyline> flowResult = fff.flow(pp, distPerVertex[i], flowIterations, finalchains[i] );
//		flowPolys.insert(flowPolys.end(), flowResult.begin(), flowResult.end());
		smoothedCurves[i] = flowResult.back();

		for (int j = 0; j < flowResult.size(); ++j)
		{
			double ratio = double(j) / flowResult.size();
			flowColors[i*flowIterations + j] = (255 * ratio, 0, 255 * (1.0 - ratio));
			flowPolys[i*flowIterations + j] = flowResult[j];
		}

        finalColors[i] = QColor(255, 0, 0);

    };
    std::cout<< "Done parallel flow \n";

    // To show the ff gradient of the polyline after the flow (replaces the flow in the ui ...)
    bool show_gradient_ff = false;

    if(show_gradient_ff)
    {
        for (MyPolyline finalCrv : smoothedCurves) {
            std::vector <MyPolyline> tmp_gradients;
            tmp_gradients = fff.get_dg_dx(finalCrv);
            for (MyPolyline e: tmp_gradients) {
                gradients.push_back(e);
            }
        }
    }

	for (int i = 0; i < flowPolys.size(); ++i)
		for (int j = 0; j < flowPolys[i].size(); ++j)
		{
			flowPolys[i][j] = Eigen::Vector2d(flowPolys[i][j].y(), flowPolys[i][j].x());
		}

    //for (int i = 191; i <= 191; ++i)
    for (int i = 0; i < smoothedCurves.size(); ++i)
		for (int j = 0; j < smoothedCurves[i].size(); ++j)
			smoothedCurves[i][j] = Eigen::Vector2d(smoothedCurves[i][j].y(), smoothedCurves[i][j].x());

    if(show_gradient_ff)
    {
        for (int i = 0; i < gradients.size(); ++i)
            for (int j = 0; j < gradients[i].size(); ++j)
                gradients[i][j] = Eigen::Vector2d(gradients[i][j].y(), gradients[i][j].x());
    }

    auto getVertexValence = [](const vertex_descriptor& v, const std::vector<std::vector<vertex_descriptor>>& allchains){
    int result = 0;
    for (auto & chain : allchains){
        result += (chain.front()==v);
        result += (chain.back()==v);
    }
    return result;
    };


//	contractDeg2(reebGraph);
    contractAllButChainEnds(reebGraph, allowedChains);
	mw.setGraph("Allowed connections", reebGraph);

    if (show_gradient_ff) mw.setPolys(gradients, flowColors,true);
    else mw.setPolys(flowPolys, flowColors, true);
//	mw.setVectorization("chopped curves", choppingResult.first);

    std::map<std::pair<int,int>, std::vector<size_t>> pixelToChainCoveringMap = buildPixelToChainsCoveringMap(reebGraph, allowedChains, origMask);
    std::vector<bool> selection;

    std::vector<MyPolyline> selectedPolys; // smooth polylines selected
    selectedPolys = selectBestPolylines(reebGraph, smoothedCurves, initialPolys, finalchains, fff, allowedChainToLabeledPointPair, pts_types, pixelToChainCoveringMap, selection);

    std::vector<std::vector<vertex_descriptor>> selectedChains;
    selectedChains.clear();
    for (size_t i_chain=0; i_chain<finalchains.size(); ++i_chain){
        if (selection[i_chain]){
            selectedChains.push_back(finalchains[i_chain]);
        }
    }

    std::vector<std::vector<double>> radii;
    std::vector<std::array<bool, 2>> protectedEnds;
    std::vector<std::array<bool, 2>> isItASpecialDeg2Vertex;
    std::vector<std::pair<PointOnCurve, PointOnCurve>> yJunctions;
    yJunctions.clear();

    for (const auto& chain : selectedChains){
        std::vector<double> tempradii;
        for (const auto& v : chain){
            tempradii.push_back(reebGraph[v].solvedWidth/2);
        }
        radii.push_back(tempradii);
//        protectedEnds.push_back({false, false});
        protectedEnds.push_back({(getVertexValence(chain.front(), selectedChains)>1), (getVertexValence(chain.back(), selectedChains)>1)});
        isItASpecialDeg2Vertex.push_back({false, false});
    }

    auto choppingResult = chopFakeEnds(selectedPolys, radii, protectedEnds, isItASpecialDeg2Vertex, yJunctions);
    std::vector<MyPolyline> finalCurves = choppingResult.first;
//    std::vector<MyPolyline> finalCurves = selectedPolys;
    std::vector<MyPolyline> allCurves;
    for (const auto& c : allowedChains)
        allCurves.push_back(vertexSeqToPolyline(c,reebGraph));
    mw.setVectorization("all polys", smoothedCurves);
    mw.setVectorization("final curves", selectedPolys);
    //mw.setPolys(finalCurves, finalColors, false);
	//mw.setPolys(allPolys);

    int last_index = finalCurves.size() -1 ;
    std::vector<Point2d> last_pts;
    for(int i = 0; i < finalCurves.size(); ++i)
    {
        for(int j = 0; j < finalCurves[i].size(); ++j)
        {
            double x = finalCurves[i][j][0];
	        double y = finalCurves[i][j][1];
	       // std::cout << i << " " << j << std::endl;
	        last_pts.push_back(Point2d (y,x));
        }
    }



    //mw.setPoints(last_pts);
    std::cout << "OK : " << last_index << std::endl;

#ifdef _WIN32
	std::cout << "All done in " << timer.elapsed()/1000.0 << " s" << std::endl;
#endif

	for (int i = 0; i < finalCurves.size(); ++i)
		finalCurves[i] = simplify(finalCurves[i], 1e-2);

	//smooth(finalCurves);
	std::string imagename = filename;
	size_t last_slash = filename.find_last_of("/\\");
	if (last_slash!=std::string::npos) {
	    imagename = filename.substr(last_slash+1);
	}

	exportSVG(svg_filename, selectedPolys, n,m, imagename);
//	exportSVG(svg_filename, smoothedCurves, n,m, imagename);
	//mw.close(); // IF YOU ONLY WANT THE SVG
	//return 0;
	if (drawMenu) {
        return a.exec();
    } else {
	    return 0;
	}
}
