#include "stdafx.h"
#include "imagequiverviewer.hpp"
#include <QtWidgets>
#include "IsLoopContractible.h"
#include "HermiteSpline.h"
ImageQuiverViewer::ImageQuiverViewer(QWidget* parent) : QLabel(parent), scale(1.0),
colors({
    { 255,102,0 },
    { 255, 146, 0 },
    { 255, 69, 0 },
    { 0, 191, 89 },
    { 209,63,240 },
    {194,42,221},
    {240,0,255},
    {84,183,44},
    {128,128,128},
}), flowHidden(true), showCircles(true), showPoints(false), polysHidden(true), showIncontractibleLoops(false) {

}

ImageQuiverViewer::~ImageQuiverViewer() {

}



void ImageQuiverViewer::setRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots)
{
	roots = newRoots;
}

void ImageQuiverViewer::setExtraRoots(const std::array<Eigen::MatrixXcd, 2>& newRoots)
{
	extraRoots = newRoots;
}

void ImageQuiverViewer::setPolys(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors, bool flow)
{
	if (!flow)
	{
		polys = newPolys;
		polyColors = colors;
	}
	else
	{
		flowPolys = newPolys;
		flowColors = colors;
	}
}

void ImageQuiverViewer::setPolys2(const std::vector<MyPolyline>& newPolys, std::vector<QColor> colors, bool flow)
{
    if (!flow)
    {
        polys = newPolys;
        polyColors = colors;
    }
    else
    {
        flowPolys = newPolys;
        flowColors = colors;
    }
}

void ImageQuiverViewer::setCovered(const std::map<vertex_descriptor, bool>& c)
{
	covered = c;
}

void ImageQuiverViewer::drawRoots(const std::array<Eigen::MatrixXcd, 2>& thoseOnes, QPainter& painter)
{
	for (int i = 0; i < thoseOnes[0].rows(); ++i)
		for (int j = 0; j < thoseOnes[0].cols(); ++j)
		{
			for (int rootIdx = 0; rootIdx < 2; ++rootIdx)
			{
				std::complex<double> vec = thoseOnes[rootIdx](i, j) * 0.5;
				QPoint origPoint((j + 0.5) * scale - vec.real() * scale, (i + 0.5) * scale - vec.imag() * scale);
				QPoint newPoint((j + 0.5) * scale + vec.real() * scale, (i + 0.5) * scale + vec.imag() * scale);
				painter.drawLine(origPoint, newPoint);
			}
		}
}

void ImageQuiverViewer::drawGraphs(QPainter& painter)
{

	for (int i = 0; i < graphs.size(); ++i)
	{
		if (!graphHidden[i])
		{
			std::map<size_t, bool> shouldDraw;
			QPen polyPen(colors[i % colors.size()]);
			polyPen.setWidth(2.0);

			QPen thickPen(colors[i % colors.size()]);
			thickPen.setWidth(5.0);

			painter.setPen(polyPen);
			QBrush brush(QColor(255, 255, 255));
			brush.setStyle(Qt::BrushStyle::SolidPattern);

			painter.setBrush(brush);
			edge_iter ei, eend;

			for (std::tie(ei, eend) = boost::edges(graphs[i]); ei != eend; ++ei)
			{
				Eigen::Vector2d edge = graphs[i][ei->m_target].location - graphs[i][ei->m_source].location;
				double len = edge.norm();
				if (false)
				{
					if (graphTags[i] != "Topo")
						//spline!
					{
						Eigen::Vector2d root1 = graphs[i][ei->m_source].root, root2 = graphs[i][ei->m_target].root;
						if (root1.norm() > 100)
							root1 = edge.normalized();
						if (root2.norm() > 100)
							root2 = edge.normalized();

						std::vector<std::array<double, 2>> possibleSigns = { { 1.0,1.0 },{ -1.0,1.0 },{ 1.0,-1.0 },{ -1.0,-1.0 } };
						double minCurv = 1e10;
						std::array<double, 2> bestSigns;
						for (auto signs : possibleSigns)
						{
							HermiteSpline hs(graphs[i][ei->m_source].location, root1 * signs[0] * len / 3, graphs[i][ei->m_target].location, root2 * signs[1] * len / 3);
							double curv = hs.integrateAbsCurvature();
							if (curv < minCurv)
							{
								bestSigns = signs;
								minCurv = curv;
							}
						}

						HermiteSpline hs(graphs[i][ei->m_source].location, root1 * bestSigns[0] * len / 3, graphs[i][ei->m_target].location, root2 * bestSigns[1] * len / 3);
						double splineLen = hs.getLength();
						std::vector<Eigen::Vector2d> pts;
						for (int i = 0; i < len * 10; ++i)
						{
							double t = i / (len * 10);
							pts.push_back(hs.getPoint(t));
						}
						for (int i = 0; i + 1 < pts.size(); ++i)
						{
							painter.drawLine(modelToScreen(pts[i]), modelToScreen(pts[i + 1]));
						}
					}
				}
				else
				{
					Eigen::Vector2d p1 = graphs[i][ei->m_source].location;
					Eigen::Vector2d p2 = graphs[i][ei->m_target].location;
					painter.drawLine(modelToScreen(p1), modelToScreen(p2));
				}
				shouldDraw[ei->m_source] = true;
				shouldDraw[ei->m_target] = true;
			}

			if (showCircles)
			{
				for (size_t v = 0; v < boost::num_vertices(graphs[i]); ++v)
				{
					if (shouldDraw[v])
					{
						if (covered[v])
							painter.setPen(thickPen);
						else
							painter.setPen(polyPen);

						Eigen::Vector2d p = graphs[i][v].location;
						painter.drawEllipse(modelToScreen(p), 3, 3);
					}
				}

				QBrush selectBrush(QColor(233, 29, 235));
				painter.setBrush(selectBrush);
				for (size_t v = 0; v < boost::num_vertices(graphs[i]); ++v)
				{
					if (covered[v])
						painter.setPen(thickPen);
					else
						painter.setPen(polyPen);

					if (graphs[i][v].clusterCurveHitSingularity && !graphs[i][v].nextToSingularity && shouldDraw[v])
					{
						Eigen::Vector2d p = graphs[i][v].location;
						painter.drawEllipse(modelToScreen(p), 5, 5);
					}
				}

				QBrush selectBrush1(QColor(55, 224, 59));
				painter.setBrush(selectBrush1);
				for (size_t v = 0; v < boost::num_vertices(graphs[i]); ++v)
				{
					if (covered[v])
						painter.setPen(thickPen);
					else
						painter.setPen(polyPen);

					if (graphs[i][v].nextToSingularity && !graphs[i][v].clusterCurveHitSingularity && shouldDraw[v])
					{
						Eigen::Vector2d p = graphs[i][v].location;
						painter.drawEllipse(modelToScreen(p), 5, 5);
					}
				}

				QBrush selectBrush2(QColor(255, 255, 59));
				painter.setBrush(selectBrush2);
				for (size_t v = 0; v < boost::num_vertices(graphs[i]); ++v)
				{
					if (covered[v])
						painter.setPen(thickPen);
					else
						painter.setPen(polyPen);

					if (graphs[i][v].nextToSingularity && graphs[i][v].clusterCurveHitSingularity && shouldDraw[v])
					{
						Eigen::Vector2d p = graphs[i][v].location;
						painter.drawEllipse(modelToScreen(p), 5, 5);
					}
				}
				painter.setBrush(brush);
			}

			if ((i == 1) && showIncontractibleLoops)
			{
				QPen polyPen(colors[(graphs.size() + 1) % colors.size()]);
				polyPen.setWidth(2.0);
				painter.setPen(polyPen);
				for (int j = 0; j < incontractibleLoops.size(); ++j)
				{
					for (int k = 0; k + 1 < incontractibleLoops[j].size(); ++k)
					{
						painter.drawLine(modelToScreen(incontractibleLoops[j][k]), modelToScreen(incontractibleLoops[j][k + 1]));
					}
				}
			}
		}
	}
}
void ImageQuiverViewer::drawPolys(QPainter& painter, const std::vector<MyPolyline>& curves, double width, std::vector<QColor> curveColors)
{
	int alpha = areGraphsVisible() ? 80 : 255;
	//QPen polyPen(QColor(0, 255, 0, alpha))
	QPen polyPen(QColor(curveColors[0].red(), curveColors[0].green(), curveColors[0].blue(), alpha));
	polyPen.setWidth(width);
	painter.setPen(polyPen);
    //painter.setBrush(QColor(curveColors[0].red(), curveColors[0].green(), curveColors[0].blue(), alpha));
    for (int i = 0; i < curves.size(); ++i)
	{
        if (curves[i].empty()) {
            continue;
        }

        if (curveColors.size() > i)
		{
            //painter.setBrush(QColor(curveColors[i].red(), curveColors[i].green(), curveColors[i].blue(), alpha));
            polyPen.setColor(curveColors[i]);
            painter.setPen(polyPen);
        }

		for (int j = 0; j < curves[i].size() - 1; ++j)
		{
			painter.drawLine(modelToScreen(curves[i][j]), modelToScreen(curves[i][j + 1]));

			//if (curves.size() < 1000000)
			//{
			//	painter.drawEllipse(modelToScreen(curves[i][j]), 3, 3);
			//}

		}
	}

	if (curves.size() < 99)
	{
        QPen polyPen(QColor(0,0,204));
        polyPen.setWidth(width);
        painter.setPen(polyPen);
		auto font = painter.font();
        font.setPointSize(16.0);
        font.setBold(true);
		painter.setFont(font);
		for (int i = 0; i < curves.size(); ++i)
		{
			if (curves[i].empty())
				continue;

			painter.drawText(modelToScreen(curves[i][curves[i].size() / 2]), QString::number(i));
		}
	}
}

void ImageQuiverViewer::drawPointSets(QPainter& painter)
{
    std::vector<std::pair<QColor, QColor>> pointStyles = {
            {QColor(0,196,196), QColor(255, 255, 255)},
            {QColor(60,255,75), {QColor(255, 255, 255)}},
            {QColor(60,255,75), {QColor(0, 0, 0)}},
            {QColor(70,240,240), {QColor(255, 255, 255)}},
            {QColor(70,240,240), {QColor(0, 0, 0)}},
            {QColor(240,50,230), {QColor(255, 255, 255)}},
            {QColor(240,50,230), {QColor(0, 0, 0)}},
            {QColor(0,120,200), {QColor(255, 255, 255)}},
            {QColor(0,0,128), {QColor(255, 255, 255)}},
            {QColor(230,25,75), {QColor(255, 255, 255)}},
    };
//    std::vector<QColor> distinctiveColors = {
//            QColor(0,196,196),
//            QColor(230,25,75),
//            QColor(70,240,240),
//            QColor(240,50,230),
//            QColor(0,0,128),
//            QColor(0,130,200),
//    };

	QPen polyPen(pointStyles.back().first);
	polyPen.setWidth(3.0);
	painter.setPen(polyPen);
	QBrush brush(pointStyles.back().second);
	brush.setStyle(Qt::BrushStyle::SolidPattern);
	painter.setBrush(brush);

	auto font = painter.font();
    font.setBold(true);
	double fontSize = std::min(10 * scale, 14.0);
	font.setPointSize(fontSize);
	painter.setFont(font);

	const double r = 4.0;

	for (int i = 0; i < pointSets.size(); ++i)
	{
//        QPen newpolyPen(distinctiveColors[i%distinctiveColors.size()]);
//        newpolyPen.setWidth(3.0);
//        painter.setPen(newpolyPen);
		if (!pointSetsHidden[i])
		{

			for (size_t v = 0; v < pointSets[i].size(); ++v)
			{
			    int pt_type = point_type[i][v];
                std::pair<QColor, QColor> style = pointStyles[pt_type%pointStyles.size()];
                QPen newpolyPen(style.first);
                newpolyPen.setWidth(3.0);
                painter.setPen(newpolyPen);
                QBrush brush(style.second);
                brush.setStyle(Qt::BrushStyle::SolidPattern);
                painter.setBrush(brush);
				Eigen::Vector2d p(pointSets[i][v].y, pointSets[i][v].x); //inverting x and y just because.
				auto pScreen = modelToScreen(p);
				painter.drawEllipse(pScreen, r, r);
				QPointF pText(pScreen.x() + 2 * r, pScreen.y() + 2 * r);
				painter.drawText(pText, QString::number(v));
			}
		}
	}
}

bool ImageQuiverViewer::areGraphsVisible()
{
	bool result = false;
	for (int i = 0; i < graphHidden.size(); ++i)
		result |= !graphHidden[i];
	return result;
}

bool ImageQuiverViewer::arePointsVisible()
{
    bool result = false;
    for (int i = 0; i < pointSets.size(); ++i)
        result |= !pointSetsHidden[i];
    return result;
}

void ImageQuiverViewer::paintEvent(QPaintEvent* event)
{
	QLabel::paintEvent(event);
	QPainter painter(this);
	int alpha = arePointsVisible() || areGraphsVisible() || !vectorizationsHidden[0] || !vectorizationsHidden[1] ? 20 + scale : 255;

	painter.setPen(QColor(216, 199, 51, alpha));
	drawRoots(roots, painter);
	painter.setPen(QColor(255, 255, 0, alpha));
	drawRoots(extraRoots, painter);

	if (!polysHidden)
	{
		double width = 1.0;
		if (polyColors.empty())
			drawPolys(painter, polys, width, { colors.back() });
		else
			drawPolys(painter, polys, width, polyColors);
	}

	drawPointSets(painter);

	if (!flowHidden)
	{
		double width = 1.0;
		if (flowColors.empty())
			drawPolys(painter, flowPolys, width, { colors.back() });
		else
			drawPolys(painter, flowPolys, width, flowColors);
	}

	drawGraphs(painter);
	drawClusters(painter);

	for (int i = 0; i < vectorizationsHidden.size(); ++i)
		if (!vectorizationsHidden[i])
			drawPolys(painter, vectorizations[i], 3.0, { colors[i] });
}

void ImageQuiverViewer::drawClusters(QPainter& painter)
{
	auto drawArrow = [&](const Eigen::Vector2d& center, const Eigen::Vector2d& dir)
	{
		auto p = painter.pen();
		QPen pen(QColor(255, 0, 0), 3);
		painter.setPen(pen);
		Eigen::Vector2d p0 = center - dir * 0.5, p1 = center + dir * 0.5;
		painter.drawLine(modelToScreen(p0), modelToScreen(p1));
		painter.drawEllipse(modelToScreen(p0), 3, 3);
	};

	QPen defaultPen(QColor(255, 255, 255));
	QPen thickPen(QColor(0, 0, 0, 100));
	thickPen.setWidth(2.0);

	if (!clustersToDraw.empty())
	{
		auto font = painter.font();
		painter.setFont(font);

		painter.setPen(defaultPen);

		for (auto cl : clustersToDraw)
		{
			QBrush br(QColor(255, 0, 0));
			painter.setBrush(br);

			for (const auto& pt : graphs[cl.first][cl.second].clusterPoints)
				painter.drawEllipse(modelToScreen(pt.p), 3, 3);

			QBrush br2(QColor(255, 255, 0, 50));
			painter.setBrush(br2);
			painter.drawEllipse(modelToScreen(graphs[cl.first][cl.second].location), graphs[cl.first][cl.second].width * scale / 2, graphs[cl.first][cl.second].width * scale / 2);

			if (graphs[cl.first][cl.second].solvedWidth > 1e-10)
			{
				QBrush br3(QColor(205, 200, 100, 50));
				painter.setBrush(br3);
				painter.setPen(thickPen);
				painter.drawEllipse(modelToScreen(graphs[cl.first][cl.second].location), graphs[cl.first][cl.second].solvedWidth * scale / 2, graphs[cl.first][cl.second].solvedWidth * scale / 2);
				painter.setPen(defaultPen);
			}

			std::set<int> adjVerts = { cl.second };
			boost::graph_traits<G>::out_edge_iterator eit, eend;
			//for (std::tie(eit, eend) = boost::out_edges(cl.second, graphs[cl.first]); eit != eend; ++eit)
			//	adjVerts.insert(eit->m_target);

			painter.setPen(colors[cl.first % colors.size()]);
			for (int v : adjVerts)
			{
				painter.drawText(modelToScreen(graphs[cl.first][v].location) + QPointF(3, 0), QString::number(v));
			}

			QPen polyPen(QColor(255, 0, 255));
			polyPen.setWidth(1.0);
			painter.setPen(polyPen);
			for (int j = 0; j + 1 < graphs[cl.first][cl.second].clusterCurve.size(); ++j)
				painter.drawLine(modelToScreen(graphs[cl.first][cl.second].clusterCurve[j]), modelToScreen(graphs[cl.first][cl.second].clusterCurve[j + 1]));
		}

		QPen pp(colors[1]);
		pp.setWidth(3.0);
		painter.setPen(pp);
		for (auto cl : clustersToDraw)
		{
			for (const auto& pt : graphs[cl.first][cl.second].clusterPoints)
			{
				if (pt.curve == graphs[cl.first][cl.second].seedCurve)
				{
					//draw root
					Eigen::Vector2d p0 = pt.p, vec = graphs[cl.first][cl.second].root;
					Eigen::Vector2d p1 = p0 + vec * 6;
					painter.drawLine(modelToScreen(p0), modelToScreen(p1));
				}
			}
		}
	}
}

void ImageQuiverViewer::mousePressEvent(QMouseEvent* event)
{
	if (event->buttons() & Qt::RightButton)
	{
		Eigen::Vector2d mousePos(event->pos().x() / scale - 0.5, event->pos().y() / scale - 0.5);
		if (event->modifiers() & Qt::ControlModifier)
		{
			std::cout << "You clicked at pixel " << (int)std::round(mousePos.y()) << " " << (int)std::round(mousePos.x()) << std::endl;
		}
		else
		{
			if (event->modifiers() & Qt::AltModifier)
				std::cout << "Selected vertices: ";
			for (int i = 0; i < graphs.size(); ++i)
			{
				if (!graphHidden[i])
				{
					for (size_t v = 0; v < boost::num_vertices(graphs[i]); ++v)
					{
						if (boost::degree(v, graphs[i]) == 0)
							continue;

						if ((graphs[i][v].location - mousePos).norm() < 3 / scale)
						{
							clustersToDraw.insert({ i,v });
							if (event->modifiers() & Qt::AltModifier)
							{
								std::cout << v << " " << " (valence " << boost::degree(v, graphs[i]) << "); ";
								oedge_iter eit, eend;
								std::cout << "Adjacent verts: ";
								for (std::tie(eit, eend) = boost::out_edges(v, graphs[i]); eit != eend; ++eit)
								{
									std::cout << eit->m_target << " ";
								}
								std::cout << std::endl;
							}
						}
					}
				}
			}
			if (event->modifiers() & Qt::AltModifier)
				std::cout << std::endl;
		}
	}
	else if (event->buttons() & Qt::LeftButton)
	{
		clustersToDraw.clear();
	}
	update();
	event->accept();
}
