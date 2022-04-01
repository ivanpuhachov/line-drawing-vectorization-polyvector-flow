#include "stdafx.h"
#include "chopThisCurve.h"
#include "intersections.h"
#include "greedyTrace.h"
#include <iostream>
std::vector<double> segmentLengths(const MyPolyline& poly)
{
	std::vector<double> result(poly.size() - 1, 0);
	for (int i = 0; i < poly.size() - 1; ++i)
		result[i] = (poly[i + 1] - poly[i]).norm();
	return result;
}

bool findFirstIntersection(const MyPolyline& poly1, const MyPolyline& poly2, double& minS)
{
	auto bbox1 = bbox(poly1);
	auto bbox2 = bbox(poly2);
	if ((bbox1[1] < bbox2[0]) || (bbox2[1] < bbox1[0]) || (bbox1[3] < bbox2[2]) || (bbox2[3] < bbox1[2]))
		return false;

	minS = std::numeric_limits<double>::max();
	for (int i=0; i<poly1.size()-1; ++i)
		for (int j = 0; j < poly2.size() - 1; ++j)
		{
			double s;
			if (get_line_intersection(poly1[i].x(), poly1[i].y(), poly1[i + 1].x(), poly1[i + 1].y(), poly2[j].x(), poly2[j].y(), poly2[j + 1].x(), poly2[j + 1].y(), nullptr, nullptr, &s, nullptr))
			{
				minS = std::min(minS, i+s);
			}
		}
	return minS < poly1.size();
}

double chopThisCurve(const std::vector<MyPolyline>& polys, int i, int whichEnd, const std::vector<std::array<bool, 2>>& protectedEnds)
{
	auto polyI = polys[i];

	if (whichEnd == 1)
		std::reverse(polyI.begin(), polyI.end());
	MyPolyline choppedEnd = polyI;

	double limit = 0;
	std::vector<double> limits(polys.size(), polyI.size());

	for (int segment = 0; segment < polys[i].size() - 1; ++segment)
	{

		for (int j = 0; j < polys.size(); ++j)
		{
			if (i == j)
				continue;
			auto b = findBoundingBox(polys[j]);
			if (b.completelyOutside(polys[i][segment], polys[i][segment + 1]))
				continue;


		}
		double s;

	}


	if (fabs(limit - polyI.size())<1e-5)
		limit = 0;

	if (whichEnd == 1)
		limit = polyI.size() - limit;

	return limit;
}
