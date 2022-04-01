#include "stdafx.h"
#include "Simplify.h"
#include <queue>

double perpendicularDistance(const Eigen::Vector2d& point, Segment segment, double& t)
{
	Eigen::Vector2d v = segment.first - segment.second;
	Eigen::Vector2d w = point - segment.second;

	Eigen::Vector2d closestPoint;

	double c1 = w.dot(v);

	if (c1<0)
	{
		closestPoint = segment.second;
		t = 1.0;
		return (closestPoint-point).norm();
	}
	double c2 = v.squaredNorm();
	if (c2 <= c1)
	{
		closestPoint = segment.first;
		t = 0.0;
		return (closestPoint - point).norm();
	}
	double b = c1 / c2;
	auto prb = segment.second + v*b;

	t = 1.0 - b;
	closestPoint = prb;
	return (prb-point).norm();
}

MyPolyline simplify(const MyPolyline & poly, double eps)
{
	if (poly.size() <= 3)
		return poly;

	std::vector<bool> mask(poly.size());
	std::queue<std::pair<int, int>> queue;
	// Find the point with the maximum distance
	queue.push({ 0,poly.size()-1});

	while (!queue.empty())
	{
		auto testPair = queue.front();
		queue.pop();

		double dmax = 0;
		int index = -1;

		for (int i = testPair.first+1; i < testPair.second; ++i)
		{
			double ignore;
			double d = perpendicularDistance(poly[i], { poly[testPair.first],  poly[testPair.second] }, ignore);
			if (d > dmax) {
				index = i;
				dmax = d;
			}
		}

		// If max distance is greater than epsilon, simplify
		if (dmax > eps) {
			queue.push({ testPair.first, index });
			queue.push({ index,testPair.second });
		}
		else {
			mask[testPair.first] = true;
			mask[testPair.second] = true;
		}

	}

	MyPolyline result;
	for (int i = 0; i < mask.size(); ++i)
	{
		if (mask[i])
			result.push_back(poly[i]);
	}
	return result;
}

MyPolyline simplify_lin_upsample(const MyPolyline & poly, double eps)
{
    if (poly.size() <= 3)
        return poly;

    int N = poly.size();
    MyPolyline new_poly;
    for(int i=1; i<N-1; ++i)
    {
        Eigen::Vector2d vec = poly[i] - poly[i-1];
        if(vec.norm() > 2 * eps)
        {
            Eigen::Vector2d uvec = vec / vec.norm();
            int segments = vec.norm() / eps;
            double seg_len = vec.norm() / segments;
            for(int j = 0; j < segments; ++j)
            {
                new_poly.push_back(poly[i-1] + j * seg_len * uvec);
            }
        }
        else
        {
            new_poly.push_back(poly[i-1]);
        }
    }
    new_poly.push_back(poly[N-1]);
    return new_poly;
}

MyPolyline simplify_lin_downsample(const MyPolyline & poly, double eps)
{
    if (poly.size() <= 3)
        return poly;

    int N = poly.size();
    MyPolyline new_poly;
    new_poly.push_back(poly[0]);
    double acc = 0;
    for(int i=1; i<N-1; ++i)
    {
        Eigen::Vector2d vec = poly[i] - poly[i-1];
        acc += vec.norm();
        if(acc > eps)
        {
            new_poly.push_back(poly[i]);
            acc = 0;
        }
    }
    new_poly.push_back(poly[N-1]);
    return new_poly;
}