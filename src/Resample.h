#include "stdafx.h"
#include "typedefs.h"

//template<class TPoint, class TVector>
MyPolyline resample(const MyPolyline & poly, int targetSize)
{
    if ((poly.size() == 1) || targetSize < 3)
        return poly;

    //std::vector<TPoint> result;
    MyPolyline result;
    //result.reserve(targetSize);

    std::vector<double> edgeLengths;
    double totalLength = 0;
    edgeLengths.reserve(poly.size());
    for (int i = 1; i < poly.size(); i++)
    {
        edgeLengths.push_back((poly[i] - poly[i - 1]).norm());
        //std::cout << "Points : " << poly[i] << " , " << poly[i-1] << "\nLength : " << edgeLengths[i-1] << std::endl;
        totalLength += edgeLengths[i - 1];
    }
    //std::cout << "Total length : " << totalLength << std::endl;

    result.push_back(poly[0]);


    for (int i = 1; i < targetSize - 1; i++)
    {
        double currentLength = 0;
        int j = 0;
        double neededLength = i*totalLength / (targetSize - 1);
        while ((j < edgeLengths.size()) && (currentLength + edgeLengths[j] < neededLength))
        {
            currentLength += edgeLengths[j];
            j++;
        }
        double alpha = (neededLength - currentLength) / edgeLengths[j];
        //TPoint p(TVector(points[j])*(1 - alpha) + TVector(points[j + 1])*alpha);
        Eigen::Vector2d p = poly[j] * (1 - alpha) + poly[j+1] * alpha;
        result.push_back(p);
    }
    result.push_back(poly.back());
    //points = result;
    return result;
}