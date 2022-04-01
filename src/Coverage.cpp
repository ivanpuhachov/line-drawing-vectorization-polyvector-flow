#include "Coverage.h"

typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> Segment;

double perpDistance(const Eigen::Vector2d& point, const Segment& segment, double& t)
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

std::set<std::pair<int, int>> computeCoverage(const G& g, const std::vector<vertex_descriptor>& path, const cv::Mat& mask, const bool includeEnds, const bool circleAround)
{
    std::set<std::pair<int, int>> result;
    for (size_t i_vertex=0; i_vertex<path.size(); ++i_vertex)
    {
        if (!includeEnds){
            if ((i_vertex==0) || (i_vertex==path.size()-1)) {
                continue;
            }
        }
        vertex_descriptor v = path[i_vertex];
        const double r = g[v].solvedWidth/2;
        if (circleAround){
            for (int i=std::floor(g[v].location.y()-r);i<std::ceil(g[v].location.y()+r); ++i)
                for (int j = std::floor(g[v].location.x() - r); j < std::ceil(g[v].location.x() + r); ++j)
                {
                    if ((mask.at<uchar>(i, j) != 0))
                        result.insert({ i,j });
                }
        }
        if (i_vertex<path.size()-1){
            vertex_descriptor v_next = path[i_vertex+1];
            const double r2 = g[v_next].solvedWidth/2;
            const double distance = (g[v].location - g[v_next].location).norm();
            if ((distance>std::max(r,r2)) || !circleAround) {
//            if ((distance>std::max(r,r2)) || true) {
                const double offset = std::max(r,r2) + distance;
                Segment s(g[v].location, g[v_next].location);
                for (int i=std::floor(g[v].location.y()-offset); i<=std::ceil(g[v].location.y()+offset); ++i)
                    for (int j = std::floor(g[v].location.x() - offset); j <= std::ceil(g[v].location.x() + offset); ++j)
                    {
                        if (i>mask.rows)
                            continue;
                        if (j>mask.cols)
                            continue;
                        if (mask.at<uchar>(i, j) != 0){
                            double t;
                            Eigen::Vector2d pt (j,i);
                            double dist = perpDistance(pt, s, t);
                            if ((t>0) and (t<1)){
                                if (dist < (1 - t) * r + t * r2) {
                                    result.insert({ i,j });
                                }
                            }
                        }
                    }
            }
        }
    }
    return result;
}

std::set<std::pair<int, int>> computeTotalCoverage(const G& g, const std::vector<std::vector<vertex_descriptor>>& existingPaths,const cv::Mat& mask)
{
    std::set<std::pair<int, int>> totalCoverage;
    for (auto path : existingPaths) {
        auto pathCoverage = computeCoverage(g,path,mask);
        totalCoverage.insert(pathCoverage.begin(), pathCoverage.end());
    }
    return totalCoverage;
}

std::set<std::pair<int, int>> computeTotalCoverageByVertices(const G& g, const std::map<vertex_descriptor,bool>& covered)
{
    std::set<std::pair<int, int>> result;
    for (size_t v=0; v<boost::num_vertices(g); ++v ){
        bool checkCoverage = false;
        try {
            checkCoverage = covered.at(v);
        } catch (std::out_of_range& e)
        {
        }
        if (checkCoverage){
            double r = g[v].solvedWidth/2;
            for (int i=std::floor(g[v].location.y()-r);i<std::ceil(g[v].location.y()+r); ++i)
                for (int j = std::floor(g[v].location.x() - r); j < std::ceil(g[v].location.x() + r); ++j)
                {
                    result.insert({ i,j });
                }
        }
    }
    return result;
}

int computeAddedCoverage(const G& g, const std::set<std::pair<int, int>>& totalCoverage, const std::vector<vertex_descriptor>& newPath, const cv::Mat& mask)
{
    std::set<std::pair<int, int>> newCoverage = computeCoverage(g, newPath, mask);
    auto common = count_if(begin(newCoverage), end(newCoverage), [&](const auto& x){ return totalCoverage.find(x) != end(totalCoverage); });
    return (newCoverage.size() - common);
}

/**
 * This function updates totalCoverage (in pixels) by adding new path on vertices
 * @param g - reebgraph
 * @param newPath - vector of edges
 * @param totalCoverage - (updated) total pixel coverage
 * @param mask - image mask
 */
void addCoverage(const G& g, const std::vector<vertex_descriptor>& newPath, std::set<std::pair<int, int>>& totalCoverage, const cv::Mat& mask)
{
    std::set<std::pair<int, int>> newCoverage = computeCoverage(g, newPath, mask);
    // DEBUG
    auto common = count_if(begin(newCoverage), end(newCoverage), [&](const auto& x){ return totalCoverage.find(x) != end(totalCoverage); });
    int newCoverageArea = newCoverage.size() - common;
    /*if (newCoverageArea < 10)
        std::cout << "WARNING: ";
    std::cout << "Adding a path, new coverage " << newCoverageArea << std::endl;*/
    // end debug
    totalCoverage.insert(newCoverage.begin(), newCoverage.end());
}

std::map<std::pair<int,int>, std::vector<size_t>> buildPixelToChainsCoveringMap(
        const G& g,
        const std::vector<std::vector<vertex_descriptor>>& existingPaths,
        const cv::Mat& mask
) {
    std::map<std::pair<int,int>, std::vector<size_t>> result;
    for (size_t i_path=0; i_path<existingPaths.size(); ++i_path) {
        std::vector<vertex_descriptor> path = existingPaths[i_path];
        std::set<std::pair<int, int>> pathCoverage = computeCoverage(g,path,mask, true, false);
        for (auto pixel_coord : pathCoverage){
            result[pixel_coord].push_back(i_path);
        }
    }
    return result;
}