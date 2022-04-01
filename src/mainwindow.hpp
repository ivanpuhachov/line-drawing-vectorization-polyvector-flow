#pragma once
#include <QWidget>
#include "ui_mainwindow.h"
#include <QtWidgets>
#include <array>
#include <Eigen/Dense>
#include "imagequiverviewer.hpp"
#include "typedefs.h"
#include "graph_typedefs.h"
#include "MyScrollArea.h"
class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();
	void setImage(QString filename, const cv::Mat& mask);
	void setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
	void setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots);
    void setPolys(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors = std::vector<QColor>(), bool flow = false);
    void setPolys2(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors = std::vector<QColor>(), bool flow = false);

    void setGraph(std::string s, const G & g);
	void setSplitVertices(const std::set<size_t>& vtx);
	void setIncontractibleLoops(const std::vector<std::vector<cv::Point2f>>& loops);
	void setVectorization(std::string s, const std::vector<MyPolyline>& vectorization);
    void setVectorization_color(std::string s, const std::vector<MyPolyline>& vectorization, std::vector<QColor> colors = std::vector<QColor>());
	void setCovered(const std::map<vertex_descriptor, bool>& covered);
	void setPointSet(std::string s, const std::vector<cv::Point2d>& pts, const std::vector<keypointType>& pts_type);
protected:
	void keyPressEvent(QKeyEvent * event);
	void redraw();
	void setScale(double newScale);
	void wheelEvent(QWheelEvent * event);
private:
	const int numGraphs, numVectorizations, numPointSets;
	Ui::MainWindow ui;
	QImage image;
	ImageQuiverViewer *imageLabel;
	MyScrollArea *scrollArea;
	double scale;
	int origWidth, origHeight;
	std::vector<QCheckBox*> checkGraphs, checkVectorizations, checkPointSets;
	QCheckBox* checkPolys, *checkIncontractibleLoops, *checkFlow;
	QRadioButton* maskRadio, *imageRadio;
	QImage actualImage, maskImage;
	std::map<std::string, int> graphNameToIdx, vectorizationNameToIdx, pointSetNameToIdx;
	QVBoxLayout *l;
};
