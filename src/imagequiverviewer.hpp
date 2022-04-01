#pragma once
#include <QtWidgets>
#include <Eigen/Dense>
#include <array>
#include "typedefs.h"
#include "graph_typedefs.h"
#include <vector>
#include "TopoGraphEmbedding.h"
class ImageQuiverViewer : public QLabel {
	Q_OBJECT

public:
	ImageQuiverViewer(QWidget * parent = Q_NULLPTR);
	~ImageQuiverViewer();

	template <typename T>
	QPointF modelToScreen(const T & p)
	{
		T result = (p + T(0.5, 0.5))*scale;
		return{ result.x(), result.y() };
	}
	void setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);;
	void setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setPolys(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors = std::vector<QColor>(), bool flow=false);
    void setPolys2(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors = std::vector<QColor>(), bool flow=false);
	void setCovered(const std::map<vertex_descriptor, bool>& covered);
	void drawRoots(const std::array<Eigen::MatrixXcd, 2>& thoseOnes, QPainter & painter);
	void drawGraphs(QPainter& painter);
	void drawPolys(QPainter & painter, const std::vector<MyPolyline>& curves, double width, std::vector<QColor> curveColors);
	void drawPointSets(QPainter& painter);
public:
	double scale;
	std::vector<bool> graphHidden, vectorizationsHidden, pointSetsHidden;
	std::vector<G> graphs;
	std::vector<std::string> graphTags;
	std::vector<std::vector<int>> graphsVertexToComponent;
	bool flowHidden, showPoints;
	bool polysHidden;
	bool showCircles;
	bool showIncontractibleLoops;
	
	std::set<size_t> splitVtx;
	std::vector<std::vector<QPointF>> incontractibleLoops;
	std::unique_ptr<Distances> tmpDistances;
	std::vector<std::vector<MyPolyline>> vectorizations;
	std::vector<std::vector<cv::Point2d>> pointSets;
	std::vector<std::vector<keypointType>> point_type;

protected:
	bool areGraphsVisible();
	bool arePointsVisible();
	void paintEvent(QPaintEvent * event);
	void drawClusters(QPainter & painter);
	virtual void mousePressEvent(QMouseEvent * event);
private:
	std::array<Eigen::MatrixXcd, 2> roots, extraRoots;
	std::vector<MyPolyline> polys, flowPolys;
	std::map<vertex_descriptor, bool> covered;
	std::set<std::pair<int, int>> clustersToDraw;
	std::vector<QColor> colors;
	std::vector<QColor> polyColors, flowColors;
};
