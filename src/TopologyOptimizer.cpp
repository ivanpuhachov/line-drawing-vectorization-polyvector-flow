#include "stdafx.h"
#include "TopologyOptimizer.h"
#include <boost/graph/filtered_graph.hpp>
#include "paal/greedy/steiner_tree_greedy.hpp"
#include "FindAllPaths.h"
#include "ShortestPathConstrained.h"
#include "Simplify.h"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <stack>

const double MINIMUM_COVERAGE_TO_BE_TAKEN_SERIOUSLY = 10;
const int MINIMUM_PIXEL_COVERAGE = 3;
const int MAXIMUM_ACTIVE_POINTS = 120;

/**
 * This function checks if the straight line between point1 and point2 lies inside the image mask.
 * Number of allowed pixels outside of the mask is thresholded by threshold.
 * @param point1
 * @param point2
 * @param origMask
 * @param threshold
 * @return
 */
bool LineLiesInMask(Eigen::Vector2d point1, Eigen::Vector2d point2, const cv::Mat& origMask, int threshold)
{
    // WN : Shouldn't we sample more ? we could miss pixels that we pass by very briefly.
    // inputs are FLIPPED (y,x) points (as in image)
    int whitepixelscounted = 0;
    double dist = (point1 - point2).norm();
    if (dist < 1) {
        return true;
    }
    int row, col;
    int prevrow = -1;
    int prevcol = -1;
    double stepx = (point1.x() - point2.x()) / std::ceil(dist);
    double stepy = (point1.y() - point2.y()) / std::ceil(dist);
    for (int i = 0; i <= std::ceil(dist); ++i) {
        if (whitepixelscounted > threshold)
            return false;

        row = int(std::round(point2.x() + i * stepx));
        col = int(std::round(point2.y() + i * stepy));
        if ((row != prevrow) || (col != prevcol)) {
            //// cv::Mat uses column-first
            if (origMask.at<uchar>(col, row) == 0) whitepixelscounted++;
        }
        prevrow = row;
        prevcol = col;
    }
    return whitepixelscounted<threshold;
}

//if once we remove these two intersecting edges, the pieces of graph we're left with contain no active points
//then we know some parts of this graph will never be traced, unless we add an intersection point
bool spuriousFrameFieldIntersection(G& g, edge_descriptor e1, edge_descriptor e2, const std::set<size_t>& activeVertices, const std::vector<cv::Point2d>& pts, double activePointDistanceThreshold)
{
    //virtually remove e1 and e2 from the graph.
    //do DFS, stopping once we found a single active vertex
    // WN : so e1 and e2 are intersecting. We return true if 2 or more of the edges endpoints cant reach any active points (within threshold)
    auto labeledPointIsClose = [&pts,&activePointDistanceThreshold](const Eigen::Vector2d& pt)

    {
        for (const auto& p : pts)
        {
            Eigen::Vector2d p_eig(p.y,p.x);
            if ((p_eig - pt).norm() < activePointDistanceThreshold)
                return true;
        }
        return false;
    };

    auto isAnyActiveVertexReachable = [&](size_t startVertex)
    {
        std::map<size_t, bool> visited;
        std::stack<size_t> q;
        q.push(startVertex);
        while (!q.empty())
        {
            size_t v = q.top();
            q.pop();

            if ((activeVertices.find(v) != activeVertices.end()) || labeledPointIsClose(g[v].location))
            {
                //active vertex is reachable!
                return true;
            }

            if (!visited[v])
            {
                visited[v] = true;
                G::out_edge_iterator e, eend;
                for (boost::tie(e, eend) = boost::out_edges(v, g);e != eend; ++e)
                {
                    if ((!visited[e->m_target]) && (*e != e1) && (*e != e2))
                    {
                        q.push(e->m_target);
                    }
                }
            }
        }
        return false;
    };

    bool e1s = isAnyActiveVertexReachable(e1.m_source), e1t = isAnyActiveVertexReachable(e1.m_target), e2s = isAnyActiveVertexReachable(e2.m_source), e2t = isAnyActiveVertexReachable(e1.m_target);

    int nNonReachable=0;
    for (auto r : { e1s, e1t, e2s, e2t })
        if (!r)
            nNonReachable++;
    return (nNonReachable >= 2);
}
/**
 * This function adds to the vertices to the graph at the intersections of the components
 * @param g
 * @param activeVertices
 * @param components
 * @param pts
 * @param activePointDistanceThreshold
 * @return
 */
std::vector<cv::Point2d> addJunctions(G& g, const std::set<size_t>& activeVertices, std::vector<int> components, const std::vector<cv::Point2d>& pts, double activePointDistanceThreshold)
{
    std::vector<cv::Point2d> result;
    auto nonAdjacent = [&g](const edge_descriptor& e1, const edge_descriptor& e2) // WN : lambda fct. [&g] is the capture clause. basically the fct has access to the reference to g.
    { // WN : checks for non adjacency and also not the same edge.
        return ((e1.m_source != e2.m_source) && (e1.m_target != e2.m_source) && (e1.m_source != e2.m_target) && (e1.m_target != e2.m_target));
    };
    // WN : For the graph g. ed : all the edges, N : nb of vertices
    auto ed = boost::edges(g);
    size_t N = boost::num_vertices(g);
    std::vector<std::pair<size_t, size_t>> edgesToAdd;
    std::vector<Eigen::Vector2d> intersectionPts;
    for (auto eit = ed.first; eit != ed.second; ++eit) // WN : unclear what's happening, is ed an iterator so ed.first is like start and ed.second is like end ? ***
    {
        for (auto eit2 = ed.first; eit2 != ed.second; ++eit2) // WN : so i think we iterate over every possible edge pair. Should verify ***
        {
            if (nonAdjacent(*eit, *eit2) && /*(components[eit->m_source] == components[eit2->m_source])*/ ((g[eit->m_source].location - g[eit2->m_source].location).norm() < 10)) //speedup. WN : check that they are non adjacent (and not the same edge too) AND not to far apart (so we can skip a TON of edge pairs)
            {
                //check intersections
                double s, t; // WN : if intersection s returns ratio of first segment to intersection point, t same for 2nd
                auto p0 = g[eit->m_source].location, p1 = g[eit->m_target].location;
                auto p2 = g[eit2->m_source].location, p3 = g[eit2->m_target].location;
                if (components[eit->m_source] == components[eit2->m_source]) {
                    if (get_line_intersection(p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), nullptr,
                                              nullptr, &s, &t)) {
                        if ((s > 1e-6) && (s < 1 - 1e-6) && (t > 1e-6) &&
                            (t < 1 - 1e-6)) // WN : I guess that means that edges intersect ?
                        {
                            if (spuriousFrameFieldIntersection(g, *eit, *eit2, activeVertices, pts,
                                                               activePointDistanceThreshold)) // WN : if the edges are too far from active points, i think ?
                            {
                                //now test that this intersection breaks the graph into for pieces, at least two of which don't have active points (i.e. will be omitted if we don't add this vertex)

                                Eigen::Vector2d intPt = p0 + (p1 - p0) * t;
                                //last check: maybe we've already added the intersection
                                bool alreadyAdded = false;
                                for (const auto &existingPt : intersectionPts) {
                                    if ((existingPt - intPt).norm() < 1e-5) {
                                        alreadyAdded = true;
                                        break;
                                    }
                                }
                                if (!alreadyAdded) // WN : queuing 4 edges to add. BUT we just created a valence = 4 vertex ????
                                {
                                    edgesToAdd.push_back({eit->m_source, N + intersectionPts.size()});
                                    edgesToAdd.push_back({eit->m_target, N + intersectionPts.size()});
                                    edgesToAdd.push_back({eit2->m_source, N + intersectionPts.size()});
                                    edgesToAdd.push_back({eit2->m_target, N + intersectionPts.size()});

                                    intersectionPts.push_back(intPt); // WN : point to add
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (auto p : intersectionPts) // WN : adding vertices to graph
    {
        size_t v = boost::add_vertex(g);
        std::cout << "Adding addJunctions vertex " << v << " at` " << p.transpose() << std::endl;
        g[v].location = p;
        g[v].clusterIdx = v;
        g[v].seedCurve = -1;
        result.push_back({ p.y(),p.x() });
    }

    for (auto e : edgesToAdd) // WN : adding edges to graph
    {
        if (!boost::edge(e.first, e.second, g).second)
        {
            auto eIdx = boost::add_edge(e.first, e.second, g);
            g[eIdx.first].weight = 0;
            //std::cout << "Adding edge " << e.first << " " << e.second << std::endl;
        }
    }

    return result; // WN : result[i] are the coord. of the i_th vertex added.
}


/**
 * This function finds the vertex in the graph, such that it is the farthest from covered vertices
 * The distance is measured by graph edges weight
 * @param g - reebGraph
 * @param coveredVertices - set of covered vertices ids
 * @param MIN_DISTANCE_THRESHOLD - distance threshold
 * @return pair (true if closer than the threshold, vertex id)
 */
std::pair<bool, size_t> findFarthestUncoveredVertexFromTheSameComponent (G& g,
                                                                         const std::set<size_t>& coveredVertices,
                                                                         const double MIN_DISTANCE_THRESHOLD
                                                                         )
{
    // add artificial node to calculate distances
    size_t artificialV = boost::add_vertex(g);
    g[artificialV].clusterIdx = artificialV;
    // add edges of weight 0 from this vertex to all other covered vertices in this component
    for (auto v: coveredVertices){
        auto newEdge = boost::add_edge(artificialV, v, g);
        g[newEdge.first].weight = 0;
    }
    // build dijkstra from artificial vertex
    std::vector<size_t> pDij(num_vertices(g)); // predecessor map for each vertex
    std::vector<double> dDij(num_vertices(g)); // distance map for each vertex
    auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    dijkstra_shortest_paths(g, artificialV,
                            predecessor_map(predMap).distance_map(distMap).weight_map(get(&Edge::weight, g)));
    // find the largest distance in dDij std::maxelem
    // mask all big values with -1 (if points are disconected, dijkstra returns 10^300
    for (size_t el=0; el<num_vertices(g); ++el){
        if (dDij[el]>100000) {
            dDij[el] = -1;
        }
    }
    auto maxelem_value_pntr = std::max_element(dDij.begin(), dDij.end());
    size_t maxelem_index = std::distance(dDij.begin(), maxelem_value_pntr);

    // remove artificial edges boost::out_edges
    boost::clear_vertex(artificialV, g);
    // remove artificial vertex
    boost::remove_vertex(artificialV, g);
    // and return the farthest point if threshold
    std::cout << "FindFarthest found vertex " << maxelem_index << ",\t at distance: " << dDij[maxelem_index]
              << ",\t at position " << g[maxelem_index].location.transpose() <<"\n";
//    return { g[maxelem_index].location.y(), g[maxelem_index].location.x() };
    return std::make_pair(dDij[maxelem_index]>MIN_DISTANCE_THRESHOLD, maxelem_index);
//        return std::make_pair((*maxelem_value_pntr>10), g[maxelem_index].location);

};

/**
 * This function finds the vertex in the graph, such that it is the farthest from covered vertices.
 * Then we calcualte pixel coverage of the path we found.
 * The distance is measured by graph edges weight
 * @param g - reebGraph
 * @param coveredVertices - set of covered vertices
 * @param currentCoverage - covered pixels coordinates
 * @param MIN_PIXEL_COVERAGE_THRESHOLD
 * @param mask
 * @return pair (bool: coverage is larger than threshold, closest vertex index)
 */
std::pair<bool, size_t> findFarthestPixelUncoveredVertexFromTheSameComponent (G& g,
                                                                              const std::set<size_t>& coveredVertices,
                                                                              const std::set<std::pair<int, int>> currentCoverage,
                                                                              const int MIN_PIXEL_COVERAGE_THRESHOLD,
                                                                              const cv::Mat& mask
)
{
    // add artificial node to calculate distances
    size_t artificialV = boost::add_vertex(g);
    // add edges of weight 0 from this vertex to all other covered vertices in this component
    for (auto v: coveredVertices){
        auto newEdge = boost::add_edge(artificialV, v, g);
        g[newEdge.first].weight = 0;
    }
    // build dijkstra from artificial vertex
    std::vector<size_t> pDij(num_vertices(g)); // predecessor map for each vertex
    std::vector<double> dDij(num_vertices(g)); // distance map for each vertex
    auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    dijkstra_shortest_paths(g, artificialV,
                            predecessor_map(predMap).distance_map(distMap).weight_map(get(&Edge::weight, g)));
    // find the largest distance in dDij std::maxelem
    // mask all big values with -1 (if points are disconected, dijkstra returns 10^300
    for (size_t el=0; el<num_vertices(g); ++el){
        if (dDij[el]>100000) {
            dDij[el] = -1;
        }
    }
    auto maxelem_value_pntr = std::max_element(dDij.begin(), dDij.end());
    size_t maxelem_index = std::distance(dDij.begin(), maxelem_value_pntr);

    std::vector<size_t> pathFound;
    size_t v = maxelem_index;
    pathFound.emplace_back(v);
    do{
        v = pDij[v];
        pathFound.emplace_back(v);
        std::cout<<v<<" ";
    } while ((pDij[v] != artificialV) && (pDij[v] != v));

    int additionalCoverage = computeAddedCoverage(g, currentCoverage, pathFound, mask);


    // remove artificial edges boost::out_edges
    boost::clear_vertex(artificialV, g);
    // remove artificial vertex
    boost::remove_vertex(artificialV, g);
    // and return the farthest point if threshold
    std::cout << "FindFarthest found vertex " << maxelem_index << ",\t at distance: " << dDij[maxelem_index]
              << ",\t at position " << g[maxelem_index].location.transpose()
              << ",\t with additional coverage " << additionalCoverage
              <<"\n";
//    return { g[maxelem_index].location.y(), g[maxelem_index].location.x() };
    return std::make_pair(additionalCoverage>MIN_PIXEL_COVERAGE_THRESHOLD, maxelem_index);
//        return std::make_pair((*maxelem_value_pntr>10), g[maxelem_index].location);

};

/**
 * Computes the distance between two vertices. Distance is measured by graph edge weight (euclid distance)
 * @param g - reebGraph
 * @param vertex1 - vertex 1 id
 * @param vertex2 - vertex 2 id
 * @return graph distance (double)
 */
double computeGraphDistanceBetweenTwoVertices(const G&g, vertex_descriptor vertex1, vertex_descriptor vertex2){
    // TODO: make it actual coverage, not euclid distance?
    std::vector<size_t> pDij(num_vertices(g)); // predecessor map for each vertex
    std::vector<double> dDij(num_vertices(g)); // distance map for each vertex
    auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx,
                                                                g)); // abstract boost container to access pDij
    dijkstra_shortest_paths(g, vertex1,
                            predecessor_map(predMap).distance_map(distMap).weight_map(get(&Edge::weight, g)));
    return dDij[vertex2];
}


/**
 * Returns the closest 2 (1 for each root) vertices of the components to the labelled points. Takes into account mask.
 * If closest 2 points are close in the graph, discard one.
 * @param p - coordinates of the labeled point to consider
 * @param component - chosen component index
 * @param reebGraph
 * @param origMask - boolean mask of the image
 * @param components - vertex_to_component vector
 * @param allRoots
 * @return array of 2 pairs of a form (closest vertex id, distance)
 */
std::array<std::pair<size_t, double>,2> findClosestVerticesInReebGraph(
        const Eigen::Vector2d& p,
        int component,
        const G& reebGraph,
        const cv::Mat& origMask,
        const std::vector<int>& components,
        const std::array<Eigen::MatrixXcd, 2>& allRoots,
        const bool doTheGraphDistanceCheck
        )
{
    // WN : Does it matter that it is a REEB graph ? This is called from a double for loop : for all components -> for all labelled points
    //returns the closest vertex _for each root_ of the frame field. WN : doesnt it have 4 roots, or root has smth to do with labelled points ?
    const int pointsOutsideTheMaskThreshold = 1;
    std::array<size_t, 2> result = { 0,0 };
    std::array<double, 2> minDist = { 1e10,1e10 };
    Eigen::Vector2d pFlip(p.y(), p.x()); // WN : Flipping x and y of the labelled pt
    std::array<int, 2> p0i = { static_cast<int>(std::round(p.x())),static_cast<int>(std::round(p.y())) }; // WN : p0i = rounded labelled pt coordinate. Is this function ever called with non integer coord pt ?
    Eigen::Vector2d root0 = _toEig(allRoots[0](p0i[0], p0i[1])).normalized(), root1 = _toEig(allRoots[1](p0i[0], p0i[1])).normalized(); // WN : the 4 roots at pixel p0i, what is the difference between root0 and root1 ?
    if (root0.norm() < 1e-5)
    {
        //warning: this means the point is outside the original mask, we resort to a single-root test
        root0 = Eigen::Vector2d(1.0, 0.0);
        root1 = Eigen::Vector2d(0, 0);
    }
    // WN: find v, the vertices (index) s.t. v is in the component received as argument
    for (size_t v = 0; v < boost::num_vertices(reebGraph); ++v)
    {
        if (components[v] != component)
            continue;

        double dist = (reebGraph[v].location - pFlip).norm(); // WN : distance between v and the labelled pt received as argument
        int idx = (abs(root0.dot(reebGraph[v].root)) > abs(root1.dot(reebGraph[v].root))) ? 0 : 1; // WN : so the index of the biggest abs( root0or1 . reebGraph[v].root). What is that last term ?
        if ((dist < minDist[idx]) && LineLiesInMask(reebGraph[v].location, pFlip, origMask, 1)) // WN : This is a bit confusing what are the TWO minDist ? it's like there's a min dist for result of above line...
            //			if (dist < minDist)
        {
            minDist[idx] = dist;
            result[idx] = v;
        }
    }
    if ((minDist[0]<1e10) && (minDist[1]<1e10)){
        if (doTheGraphDistanceCheck){
            double dist = computeGraphDistanceBetweenTwoVertices(reebGraph, result[0], result[1]);
            std::cout << "Those two points " << result[0] << " " << result[1] <<  " are this close in graph (" <<  dist << ") \n";
            if (dist < 1.5 * MINIMUM_COVERAGE_TO_BE_TAKEN_SERIOUSLY) {
                // if two artificial points are too close in the graph, pick 1 and discard the other
                std::cout << "Discarding the other one, that is the farthest from the labeled point ";
                int other_index = minDist[0] < minDist[1] ? 1 : 0;
                std::cout << result[other_index] << "\n";
                minDist[other_index] = 100;
            }
        }
    }
    //        std::cout << "Point (" << p.transpose() << ") line to component "<< component << " has " << std::endl << LineLiesInMask(reebGraph[result].location, pFlip, true) << std::endl;
    return { std::make_pair(result[0], minDist[0]), std::make_pair(result[1], minDist[1]) };
};

/**
 * This function computes connected components.
 * See https://www.boost.org/doc/libs/1_68_0/libs/graph/doc/connected_components.html
 * boost::connected_components() functions compute the connected components of an undirected graph using a DFS-based approach.
 * A connected component of an undirected graph is a set of vertices that are all reachable from each other.
 * @param reebGraph
 * @return pair of (vertex_to_component vector, component_to_vertices_vector vector)
 */
std::pair< std::vector<int>, std::vector<std::vector<vertex_descriptor>>> computeComponents(G& reebGraph)
{
    std::vector<int> vertex_to_component(boost::num_vertices(reebGraph)); // WN : components is a vector of lenght num_vertices(reeb_Graph)
    int nComponents = boost::connected_components(reebGraph, &vertex_to_component[0]); // WN : nComponents is the number of components, in the reeb graph, it's the number of connected "subgraph". This line also fills in the components vector. components[v] is the component nb. of a given vertex v of the graph.
    std::vector<std::vector<vertex_descriptor>> component_to_vertices_vector(nComponents); // WN : Kind of the inverse of the components vector. myComponents[i] contains a vector of all the vertices of component i

    for (size_t v = 0; v < boost::num_vertices(reebGraph); ++v)
        component_to_vertices_vector[vertex_to_component[v]].push_back(v);
    return { vertex_to_component, component_to_vertices_vector};
}

void moveGraphVertex(
        G& reebGraph,
        const cv::Point2d& point,
        const size_t vertexToMove
        ) {
    // here we move the vertex to labeled point location
    Cluster c = reebGraph[vertexToMove];
    c.location = Eigen::Vector2d(point.y, point.x);
    c.split = false;
    c.clusterIdx = reebGraph[vertexToMove].clusterIdx;
    c.width = reebGraph[vertexToMove].width;
    c.solvedWidth = reebGraph[vertexToMove].solvedWidth;
    reebGraph[vertexToMove] = c;
}


/**
 * Compute a vector from component id to a set of (labeled point, corresponding closest point in the component).
 * Closest points are thresholded by internal activePointDistanceThreshold.
 * Closest points are computed w.r.t. frame field directions.
 * if (moveVertices) changes the coordinates of the closest points to the exact position of the labeled point.
 * @param reebGraph (may be changed)
 * @param pts - vector coordinates of labeled points
 * @param origMask - mask
 * @param nComponents - number of connected components in the graph
 * @param components - vertex_to_component vector
 * @param allRoots - roots of the framefield
 * @param moveVertices - boolean flag to change vertex coordinates
 * @param doTheGraphDistanceCheck - boolean flag to check the distance between active vertices
 * @return vector from component id to a set of (labeled point, corresponding closest point in the component)
 */
std::vector<std::set<std::pair<int, int>>> computeActivePointsPerComponent(
        G& reebGraph,
        const std::vector<cv::Point2d>& pts,
        const cv::Mat& origMask,
        int nComponents,
        const std::vector<int>& components,
        const std::array<Eigen::MatrixXcd, 2>& allRoots,
        const bool moveVertices,
        const bool doTheGraphDistanceCheck
        )
{
    // WN : For each component, call findClosestVerticesInReebGraph (which returns pairs) for every labeled pt then iterate over those pairs (still for each component) : if the second element is smaller than threshold : do activePtsPerComponent[component].insert({ i, it.first }); (activePtsperComponent[component] is a set of pairs of which the first element is the labelled point index and the second is the
    const double activePointDistanceThreshold = 4; // max distance in pixel between labeled point and graph vertex
    //stores (idx of point, idx of closestVertex)
    std::vector<std::set<std::pair<int, int>>> activePtsPerComponent(nComponents); // WN : activePtsPerComponent[a] = {{b, c}, ...} where a = component, b = labelled pt, c = closest point of component (if close enough)

    for (int component = 0; component < nComponents; ++component)
    {

        for (int i = 0; i < pts.size(); ++i) {
            auto vStartPairs = findClosestVerticesInReebGraph(Eigen::Vector2d(pts[i].x, pts[i].y), component, reebGraph, origMask, components, allRoots, doTheGraphDistanceCheck);
            for (auto it : vStartPairs) // WN : iterate over the 2 closest vertices returned by above fct. they are pairs. it.first is the vertex index, it.second is the distance
            {
                if (it.second < activePointDistanceThreshold) { // if distance between labeled point and graph vertex is less than a threshold
                    //                std::cout << "Point (" << i << ") line to component "<< component << " has " << std::endl << LineLiesInMask(reebGraph[vStartPair.first].location, Eigen::Vector2d(pts[i].y, pts[i].x), true) << std::endl;
                    //std::cout << "Active vertex: " << it.first << std::endl;
                    activePtsPerComponent[component].insert({ i, it.first });

                    if (moveVertices)
                    {
                        moveGraphVertex(reebGraph, pts[i], it.first);
                    }
                }
            }
        }
    }
    return activePtsPerComponent;
}


std::set<std::pair<int, int>> filterActivePtsByDistance(
        G& reebGraph,
        std::set<std::pair<int, int>> currentActivePts, // vector of (b,c)  b = labelled pt, c = closest vertex to b
        const std::vector<cv::Point2d>& pts, // labeled points positions
        const bool moveVertices // move the vertices location to labeled point location
        ) {
    std::set<std::pair<int, int>> result;

    // handy function to insert active point pair (b,c)  b = labelled pt, c = closest vertex to b
    auto addPairToResult = [&reebGraph, &pts, &moveVertices, &result](std::pair<int, int> pairtoinsert)
    {
        result.insert(pairtoinsert);
        if (moveVertices) {
            moveGraphVertex(reebGraph, pts[pairtoinsert.first], pairtoinsert.second);
        }
    };

    while (!currentActivePts.empty()) {
        std::set<std::pair<int, int>>::iterator elem = currentActivePts.begin();
        std::pair<int, int> beginpair = std::make_pair(elem->first, elem->second);
        currentActivePts.erase(elem); // remove elem from set
        if (!currentActivePts.empty()){
            // find second active vertex associated with this labeled point
            std::set<std::pair<int, int>>::iterator elem2 = std::find_if(
                    currentActivePts.begin(),
                    currentActivePts.end(),
                    [&elem](const std::pair<int, int> el) {return (elem->first==el.first); }
                    );
            if (elem2!=currentActivePts.end()) {
                //if we find second vertex associated with labeled point, compute distance between them
                double dist = computeGraphDistanceBetweenTwoVertices(reebGraph, elem->second, elem2->second);
                std::cout << "Those two points " << elem->second << " " << elem2->second <<  " are this close in graph (" <<  dist << ") \n";
                std::pair<int, int> secondpair = std::make_pair(elem2->first, elem2->second);
                if (dist < 1.5 * MINIMUM_COVERAGE_TO_BE_TAKEN_SERIOUSLY) {
                    std::cout << "Discarding the one, that is the farthest from the labeled point ";
                    // computing distance between vertex and labeled point location
                    Eigen::Vector2d pFlip(pts[elem->first].y, pts[elem->first].x);
                    double dist = (reebGraph[elem->second].location - pFlip).norm();
                    double dist2 = (reebGraph[elem2->second].location - pFlip).norm();
                    if (dist < dist2) {
                        addPairToResult(beginpair);
                    } else {
                        addPairToResult(secondpair);
                    }
                } else {
                    std::cout << "They are sufficiently far away\n";
                    addPairToResult(beginpair);
                    addPairToResult(secondpair);
                }
            } else {
                // if 2nd vertex not found, just include elem in the result
                addPairToResult(beginpair);
            }
        } else {
            // if that was the last point in set, include it in the result
            addPairToResult(beginpair);
        }
    }
    return result;
}



std::vector<MyPolyline> splitPolyAtActivePts(MyPolyline poly, std::set<std::pair<int, int>> activePts, G& reebGraph)
{
    std::vector<MyPolyline> polys;
    bool doNothing = false;
    if(!doNothing)
    {
        MyPolyline tmp_poly;
        bool foundAHit = false;
        for (auto &p : poly) {
            tmp_poly.push_back(p);
            if (p != poly.back() && p != poly.front()) {
                for (auto &actPt : activePts) {
                    auto& loc = reebGraph[actPt.second].location;
                    if (p == loc && actPt.first != -1) {
                        std::cout << "found a hit ! ! " << std::endl;
                        foundAHit = true;
                        break;
                    }
                }
            }
            if (foundAHit) {
                polys.push_back(tmp_poly);
                tmp_poly.clear();
                tmp_poly.push_back(p);
                foundAHit = false;
            } else if (p == poly.back())
                polys.push_back(tmp_poly);
        }
        return polys;
    }
}

void untieHighValenceVertex(std::vector<std::vector<edge_descriptor>>& chains, std::vector<MyPolyline>& polys, std::vector<size_t> adjacentChainsIdx, G& g)
{
    std::cout<<"Inside untie fct"<<std::endl;

    // Here's the plan. We will transform a valence n vertex + it's n chains into n - 1 chain, of which n - 2 will be stems.
    // First idea : the longest chain will not be a stem. This chain will be replaced by a merge of itself and the smallest edge (because it means shortest graph distance between between chain ends)

    auto getEdgeLength = [&g](const edge_descriptor e)
    {
        double l = (g[e.m_target].location - g[e.m_source].location).norm();
        return l;
    };

    for (int i = 0 ; i<adjacentChainsIdx.size();i++)
    {
        std::cout<<"chain "<<adjacentChainsIdx[i]<<std::endl;
        printChain(chains[adjacentChainsIdx[i]], g, true);
    }

    std::vector<double> lengthChains;
    for (int cIdx : adjacentChainsIdx)
    {
        double len = 0.0;
        for (auto &e : chains[cIdx])
        {
            len += getEdgeLength(e);
        }
        lengthChains.push_back(len);
    }

    std::map<size_t, std::vector<edge_descriptor>> newChain;
    int minElementIndex;
    int minAdjChainIdx;
    int maxElementIndex = std::max_element(lengthChains.begin(),lengthChains.end()) - lengthChains.begin();
    int maxAdjChainIdx = adjacentChainsIdx[maxElementIndex];
    lengthChains.erase(lengthChains.begin() + maxElementIndex);
    adjacentChainsIdx.erase(adjacentChainsIdx.begin() + maxElementIndex);

    // Since I dont get how to iterate over values of map in order of sorted values, have to do some weird parallel vector management with .erase() ....
    // Merges 2 chains at a time, keeps them for assignation into the chains variable for after that for loop.
    int nbLoop = adjacentChainsIdx.size();
    for (int i = 0; i<nbLoop; i++)
    {
        minElementIndex = std::min_element(lengthChains.begin(),lengthChains.end()) - lengthChains.begin();
        minAdjChainIdx = adjacentChainsIdx[minElementIndex];
        lengthChains.erase(lengthChains.begin() + minElementIndex);
        adjacentChainsIdx.erase(adjacentChainsIdx.begin() + minElementIndex);

        std::cout<<"Merging chain "<<maxAdjChainIdx<<" and "<<minAdjChainIdx<<std::endl;
        newChain[maxAdjChainIdx] = mergeChains(chains[maxAdjChainIdx], chains[minAdjChainIdx], g);

        maxAdjChainIdx = minAdjChainIdx;
    }

    // assignation of values computed in previous loop
    for (auto& cPlusP : newChain)
    {
        chains[cPlusP.first] = cPlusP.second;
    }

    // since we promised n-1 chains, that's the one we will delete.
    chains.erase(chains.begin() + maxAdjChainIdx);
    polys.erase(polys.begin() + maxAdjChainIdx);
}

/**
 * This function iterates through chains and finds chain endpoints that are not active and are not close to the labeled point.
 * Distance is thresholded by internal activePointDistanceThreshold.
 * @param g - reebGraph
 * @param labeledpoints - coordinates of the labeled point
 * @param chains
 * @param activePtsPerComponent = vector to
 * @param selectedComponentId
 * @param activePointDistanceThreshold
 * @return
 */
std::map<size_t, std::set<size_t>> findUnexpectedEndChainsVertices(
        const G& g,
        const std::vector<cv::Point2d>& labeledpoints,
        const std::vector<std::vector<edge_descriptor>>& chains,
        std::vector<std::set<std::pair<int, int>>> activePtsPerComponent,
        const int selectedComponentId,
        const double activePointDistanceThreshold=0.0)
        {
    std::map<size_t, std::set<size_t>> nonValence2UnepectedVertexToAdjacentChainIds;
    for (int cIdx=0; cIdx<chains.size(); cIdx++) // WN : iterate over chains. cIdx is current chain idx, c is current chain
    {
        const auto &c = chains[cIdx];
        for (auto v : { c.front().m_source, c.back().m_target })  // WN : iterate over beginning and ending of chain. v is that vertex
        {
            //check if current vertex is an active vertex for this component
            if (std::find_if(activePtsPerComponent[selectedComponentId].begin(),
                             activePtsPerComponent[selectedComponentId].end(),
                             [&v](const std::pair<int, int>& el) {return el.second == v; }) == activePtsPerComponent[selectedComponentId].end())
            {
                bool foundAnnotationPt = false;
                for (int i = 0; i < labeledpoints.size(); ++i) // WN : Iterate over all labelled pts.
                {
                    Eigen::Vector2d pFlip(labeledpoints[i].y, labeledpoints[i].x);
                    if ((pFlip - g[v].location).norm() < activePointDistanceThreshold) // WN : If the labelled pt is close enough to this end of the chain. We found an Annotation pt.
                    {
                        foundAnnotationPt = true;
                        std::cout << "Found annotation point: " << v << " at location " << g[v].location.transpose() << std::endl;
                        break;
                    }
                }

                if (!foundAnnotationPt)
                {
                    nonValence2UnepectedVertexToAdjacentChainIds[v].insert(cIdx); // WN : If we didn't find an annotation pt, it means the this end of the chain might have a valence > 2. So we start accumulating the chain indices.
                }
            }
        }
    }
    return nonValence2UnepectedVertexToAdjacentChainIds;
}


bool checkChainOfChains (
        const std::vector<std::pair<std::pair<size_t, size_t>, std::vector<edge_descriptor>>> in,
        const std::vector<std::pair<std::pair<size_t, size_t>, std::vector<edge_descriptor>>> notIn,
        const std::set<size_t> allAdjChains,
        std::vector<std::vector<edge_descriptor>>& finalChoice,
        const G& g)
{
    std::cout<<std::endl<<"Begin checkChainOfChains ";
//    std::cout<<"in chains (2 int each)";
//    for (auto e : in)
//    {
//        std::cout<<e.first.first<<" "<<e.first.second<<" - ";
//    }
//    std::cout<<"not in chains (2 int each)";
//    for (auto e : notIn)
//    {
//        std::cout<<e.first.first<<" "<<e.first.second<<" - ";
//    }
//    std::cout<<std::endl;
    if (in.size() == 0)
    {
        std::cout << " in size = 0" << std::endl;
        for (int i = 0; i < notIn.size(); i++)
        {
            auto newIn = in;
            auto newNotIn = notIn;
            newIn.push_back(notIn[i]);
            newNotIn.erase(newNotIn.begin() + i);
            if (checkChainOfChains(newIn, newNotIn, allAdjChains, finalChoice, g)) return true;
        }
        std::cout<<"Didnt find any OK combination"<<std::endl;
        return false;
    }
    //otherwise, test for intersection first
    for (auto e1 : in)
    {
        auto c1 = e1.second;
        auto p0 = g[c1.front().m_source].location;
        auto p1 = g[c1.back().m_target].location;
        for (auto e2 : in)
        {
            if (e1 == e2) continue;
            auto c2 = e2.second;
            auto p2 = g[c2.front().m_source].location;
            auto p3 = g[c2.back().m_target].location;
            double s, t;
            if (get_line_intersection(p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), nullptr ,nullptr ,&s,&t)) {
                if ((s > 1e-6) && (s < 1 - 1e-6) && (t > 1e-6) && (t < 1 - 1e-6)) {
                    std::cout << "Intersecting" << std::endl;
                    return false;
                }
            }
        }
    }
    std::cout<<"No intersection found"<<std::endl;

    std::vector<std::pair<size_t, size_t>> basicChains;
    size_t begin, end;
    std::set<size_t> needToCover = allAdjChains;

    for (auto e : in) basicChains.push_back(e.first);
    begin = basicChains[0].first;
    end = basicChains[0].second;
    basicChains.erase(basicChains.begin());
    needToCover.erase(begin);
    needToCover.erase(end);

    bool stillPossible = true;

    // Recap : basicChains empty first : not enough chains -> need more chains,
    // needToCover empty first : we cover all (very good),
    // stillPossible = false : the big chain we are (conceptually) building is stuck (think of islands !) -> need more chains
    // begin == end : we have a cycle (no good)
    while (basicChains.size() != 0 && needToCover.size() != 0 && stillPossible)
    {
        if (begin == end) // a loop
        {
            stillPossible = false;
            break;
        }
        bool foundMatch = false;
        for (int i = 0; i < basicChains.size(); i++)
        {
            if (basicChains[i].first == begin && needToCover.find(basicChains[i].second) != needToCover.end()) // if we continue our conceptual big chain AND we dont go back to an already covered vertex
            {
                foundMatch = true;
                begin = basicChains[i].second;
                basicChains.erase(basicChains.begin() + i);
                needToCover.erase(basicChains[i].second);
                break;
            }
            else if(basicChains[i].first == end && needToCover.find(basicChains[i].second) != needToCover.end())
            {
                foundMatch = true;
                end = basicChains[i].second;
                basicChains.erase(basicChains.begin() + i);
                needToCover.erase(basicChains[i].second);
                break;
            }
            else if (basicChains[i].second == begin && needToCover.find(basicChains[i].first) != needToCover.end())
            {
                foundMatch = true;
                begin = basicChains[i].first;
                basicChains.erase(basicChains.begin() + i);
                needToCover.erase(basicChains[i].first);
                break;
            }
            else if(basicChains[i].second == end && needToCover.find(basicChains[i].first) != needToCover.end())
            {
                foundMatch = true;
                end = basicChains[i].first;
                basicChains.erase(basicChains.begin() + i);
                needToCover.erase(basicChains[i].first);
                break;
            }
        }
        if (!foundMatch) stillPossible = false;
    }

    if (begin == end)
    {
        std::cout<<"There's a loop, these chains wont work"<<std::endl;
        return false;
    }

    else if (needToCover.size() == 0)
    {
        std::cout<<"Found a good combination"<<std::endl;
        for (auto e : in) {
            std::cout << e.first.first << " " << e.first.second << std::endl;
            finalChoice.push_back(e.second);
        }
        return true;
    }

    else if(!stillPossible || basicChains.size() == 0) {
        if (notIn.size() == 0) {
            std::cout << "No more chains to add, dead end" << std::endl;
            return false;
        }
        else {
            std::cout<< "Not enough chains, lets add some"<<std::endl;
            for (int i = 0; i < notIn.size(); i++) {
                auto newIn = in;
                auto newNotIn = notIn;
                newIn.push_back(notIn[i]);
                newNotIn.erase(newNotIn.begin() + i);
                if (checkChainOfChains(newIn, newNotIn, allAdjChains, finalChoice, g)) return true;
            }
            std::cout<<"In this solution branch, no OK combination found"<<std::endl;
            return false;
        }
    }


}
/**
 * This function finds the labeled point id connected to selected point in the component
 * @param activePtsPerComponent
 * @param selectedPointIdx
 * @param selectedComponentIdx
 * @return labeled point id, -1 if not found
 */
int getTheLabeledPointConnectedToMe (
        const std::vector<std::set<std::pair<int, int>>> &activePtsPerComponent,
        const size_t selectedPointIdx,
        const int selectedComponentIdx
        ){
    // if labeled point not found, then -1
    int foundLabelIdx = -1;
    for (auto pair : activePtsPerComponent[selectedComponentIdx]){
        if (pair.second == selectedPointIdx){
            foundLabelIdx = pair.first;
            break;
        }
    }
    return foundLabelIdx;
};

/**
 * This function finds the labeled point id connected to selected point and looks in all components
 * @param activePtsPerComponent - vector component to the set of (labeled point id, vertex id)
 * @param selectedPointIdx
 * @return
 */
int getLabeledPointConnected (
        const std::vector<std::set<std::pair<int, int>>> &activePtsPerComponent,
        const size_t selectedPointIdx
        ){
    // if labeled point not found, then -1
    int foundLabelIdx = -1;
    for (int i_component=0; i_component<activePtsPerComponent.size(); i_component++) {
        for (auto pair : activePtsPerComponent[i_component]) {
            if (pair.second == selectedPointIdx) {
                foundLabelIdx = pair.first;
                break;
            }
        }
    }
    return foundLabelIdx;
};

/**
 * This function builds a steiner tree on activepoints and extracts chains
 * @param reebGraph
 * @param activePtsPerComponent
 * @param selectedComponentId
 * @param activepoints
 * @return
 */
std::vector<std::vector<edge_descriptor>> buildChainsFromSteinerTree(
        const G& reebGraph,
        const std::set<size_t>& activepoints
        ){
    std::set<edge_descriptor> steinerEdges; // WN : preparing to build a steiner tree.
    std::vector<int> color(num_vertices(reebGraph)); // WN : preparing a graph coloring.

    {
        auto c = &color[0];
        for (auto v : activepoints)
            boost::put(c, v, paal::Terminals::TERMINAL); // WN : Not sure .. We are coloring the active pts differently. Are they all tagged as terminals ?
    }
    // WN : Steiner tree building. TODO Read doc
    auto index = get(boost::vertex_index, reebGraph);
    paal::steiner_tree_greedy(
            reebGraph, std::inserter(steinerEdges, steinerEdges.begin()),
            boost::vertex_color_map(
                    boost::make_iterator_property_map(color.begin(), index)).weight_map(boost::get(&Edge::weight, reebGraph)));

    struct EdgeInSteinerTree // WN : Unsure what is the use of this struct, the method I guess ?
    {
        EdgeInSteinerTree() {}
        EdgeInSteinerTree(const std::set<edge_descriptor>& steinerEdges) : steinerEdges(steinerEdges) {}
        bool operator()(const edge_descriptor& e) const { // WN : This is to check if and edge e is part of the steiner tree.
            return (steinerEdges.find(e) != steinerEdges.end());
        }
        std::set<edge_descriptor> steinerEdges;
    };

    EdgeInSteinerTree filter(steinerEdges);
    boost::filtered_graph<G, EdgeInSteinerTree> fg(reebGraph, filter); // WN : fg is a filtered graph TODO read doc
    std::map<edge_descriptor, size_t> myChain;
    std::vector<std::vector<edge_descriptor>> chainsFromSteiner = chainDecompositionWithActive(fg, activepoints, myChain); // chains contain edge chains in this component
    return chainsFromSteiner;
}

/**
 * This function updates map, marking some vertices as covered if they lie close to provided chain.
 * Distance is computed as perpendicular distance to the edge and is compared to the coverage radiuses of vertices.
 * @param reebGraph
 * @param chainToProcess
 * @param componentToVertices
 * @param selectedComponentId
 * @param vertexIsCoveredByRadius
 */
void markVerticesCoveredByRadius (
        const G& reebGraph,
        const std::vector<edge_descriptor> &chainToProcess,
        const std::vector<std::vector<vertex_descriptor>> &componentToVertices,
        const int selectedComponentId,
        std::map<vertex_descriptor, bool>& vertexIsCoveredByRadius
        ){
    for (const edge_descriptor &e : chainToProcess)
    {
//                        vertexIsCoveredByRadius[e.m_source]=true;
        Segment s(reebGraph[e.m_source].location,
                  reebGraph[e.m_target].location);
        double r1 = 0.5 * reebGraph[e.m_source].solvedWidth, r2 =
                0.5 * reebGraph[e.m_target].solvedWidth;
        for (size_t v : componentToVertices[selectedComponentId])
        {
            double t;
            auto pt = reebGraph[v].location;
            double dist = perpendicularDistance(pt, s, t);
            if (dist < (1 - t) * r1 + t * r2)
                vertexIsCoveredByRadius[v] = true;
        }
    }
//                    vertexIsCoveredByRadius[chain.back().m_target]=true;
}


std::set<size_t> getVerticesUsedInChains (
        const std::vector<std::vector<edge_descriptor>>& chains
        ) {
    std::set<size_t> result;
    std::cout << "Getting covered vertices from " << chains.size() << " chains \n";
    for (const auto &chain : chains) {
        result.insert(chain.front().m_source);
        for (const auto &e : chain) {
            result.insert(e.m_target);
        }
    }
    return result;
};


void tryToAddNewChain (
        const G& reebGraph,
        const cv::Mat& origMask,
        const std::vector<edge_descriptor>& newChain,
        std::vector<std::vector<edge_descriptor>>& chains,
        std::set<std::pair<int, int>>& TotalPixelCoverage
        ){
    std::cout << "We want to add path ";
    printChain(newChain, reebGraph, true);

    // find what new edges are we really adding (subchains)
    std::vector<edge_descriptor> subchain;
    std::vector<std::vector<edge_descriptor>> newsubchains;
    bool newSubchainIsBuilding = false;
    // extract only new subchains from newChain
    for (auto candidateEdge : newChain){
        // check if this edge is already in chains
        if (!(checkIfEdgeIsInChains(candidateEdge, chains))) {
            if (newSubchainIsBuilding) {
                // add to building subchain
                subchain.push_back(candidateEdge);
            } else {
                // start building a new subchain
                newSubchainIsBuilding = true;
                subchain.push_back(candidateEdge);
            }
        } else {
            if (newSubchainIsBuilding) {
                // finish building the subchain
                newSubchainIsBuilding = false;
                newsubchains.push_back(subchain);
                subchain.clear();
            }
        }
    }
    // if the last edge in newChain is uncovered too - finish building the subchain
    if (newSubchainIsBuilding) {
        newSubchainIsBuilding = false;
        newsubchains.push_back(subchain);
        subchain.clear();
    }
    std::cout << "I found this uncovered subchains: \n";
    for (auto & c : newsubchains) {
        printChain(c, reebGraph);
    }

    // select endpoints of new subchains
    std::set<vertex_descriptor> newendpoints;
    for (auto & subchain : newsubchains) {
        newendpoints.insert(subchain.front().m_source);
        newendpoints.insert(subchain.back().m_target);
    }

    // add subchains to the coverage
    for (auto & subchain : newsubchains) {
        chains.push_back(subchain);
        addCoverage(reebGraph, chainToVerticesSeq(subchain), TotalPixelCoverage, origMask);
    }

    // split existing chains by selected endpoints
    for (const vertex_descriptor & splitvertex : newendpoints){
        for (size_t i_chain=0; i_chain<chains.size(); ++i_chain) {
            if (vertexIsInsideChain(splitvertex, chains[i_chain])){
                auto splittingresult = splitChainByVertex(splitvertex, chains[i_chain]);
                chains[i_chain] = splittingresult.first;
                chains.push_back(splittingresult.second);
            }
        }
    }

    std::cout << "Now chains are " << chains.size() << "\n";
    for (int i = 0; i<chains.size(); i++)
    {
        std::cout << "chain : " << i;
        printChain(chains[i], reebGraph, true);
    }

    // FOLLOWING FOR LOOP : IS FOR FIXING Y JUNCTIONS AFTER ADDING A COVERING POLYLINE
    bool fixYJunctionsNewEdge = false;
    if (fixYJunctionsNewEdge)
    {
        for (size_t cIdx = 0; cIdx < chains.size(); cIdx++)
        {

            int yJunctionIdx = -1;

            std::cout << "Step 1 : Looking for val 3 vertices in newChain and chain " << cIdx
                      << std::endl;
            std::vector<int> valences(boost::num_vertices(reebGraph),
                                      0); // index is vertex index, value is valence of vertex
            std::set<edge_descriptor> uniqueEdges;

            //Printing the chain and also building uniqueEdges
            std::cout << "Front of chain : " << chains[cIdx].front().m_source << " : "
                      << reebGraph[chains[cIdx].front().m_source].location.x() << ", "
                      << reebGraph[chains[cIdx].front().m_source].location.y() << std::endl;

            for (auto &e : chains[cIdx]) {
                uniqueEdges.insert(e);
                //    std::cout<<e.m_source<<std::endl;
            }
            std::cout << "back of chain : " << chains[cIdx].back().m_target << " : "
                      << reebGraph[chains[cIdx].back().m_target].location.x() << ", "
                      << reebGraph[chains[cIdx].back().m_target].location.y() << std::endl;
            for (auto &e : newChain) {
                uniqueEdges.insert(e);
            }
            for (auto &e : uniqueEdges) {
                valences[e.m_source]++;
                valences[e.m_target]++;
            }

            bool alreadyFoundJunction = false;
            for (int i = 0; i < valences.size(); i++)
                if (valences[i] == 3) {
                    std::cout << "   vertex " << i << " has valence : " << valences[i]
                              << ". Coords : " << reebGraph[i].location.x() << ", "
                              << reebGraph[i].location.y() << std::endl;
                    if (alreadyFoundJunction) yJunctionIdx = -2;
                    else if (i == chains[cIdx].front().m_source ||
                             i == chains[cIdx].back().m_target || i == newChain.front().m_source ||
                             i == newChain.back().m_target)
                        std::cout << "hit with end pt" << std::endl;
                    else {
                        yJunctionIdx = i;
                        alreadyFoundJunction = true;
                    }
                }

            std::cout << "end of step 1. Chain " << cIdx << " with y junction idx : " << std::endl;
            if (yJunctionIdx == -1) {
                std::cout << "Didnt find a y junction" << std::endl;
                continue;
            } else if (yJunctionIdx == -2) {
                std::cout << "Found 2 y junctions for one chain." << std::endl;
                continue;
            } else std::cout << yJunctionIdx << std::endl<<std::endl;

            std::cout<<"Printing chain : "<<cIdx<<std::endl;

            //printChain(chains[cIdx], reebGraph);
            //printChainByEdges(chains[cIdx], reebGraph);



            std::cout << "Step 2 : 2 chains -> 3 chains meeting at y junction";


            std::vector<edge_descriptor> c1;
            std::vector<edge_descriptor> c2;
            std::vector<edge_descriptor> c3;
            std::vector<edge_descriptor> c4;
            bool buildingC1 = true;
            for (auto &e : chains[cIdx]) {
                if (e.m_source == yJunctionIdx) {
                    //std::cout << "origin pt reached" << std::endl;
                    buildingC1 = false;
                }
                if (buildingC1) {
                    c1.push_back(e);
                } else {
                    c2.push_back(e);
                }
            }

            bool buildingC3 = true;
            for (auto &e : newChain) {
                if (e.m_source == yJunctionIdx) {
                    buildingC3 = false;
                }
                if (buildingC3) {
                    c3.push_back(e);
                } else {
                    c4.push_back(e);
                }
            }
            // Now we have to choose between c3 and c4. If c3 shares and edge with c1 or c2, choose c4.
            bool c3Good = true;
            for (auto &e : c3) {
                for (auto &e1 : c1) {
                    if (e == e1) {
                        c3Good = false;
                    }
                }
                for (auto &e2 : c2) {
                    if (e == e2) {
                        c3Good = false;
                    }
                }
            }
            if (!c3Good) c3 = c4;

            std::cout << "End of step 2 : building the 3 chains." << std::endl << std::endl;
            std::cout << "    chain1 :" << std::endl;
            std::cout << "index : " << c1.front().m_source << " coords : "
                      << reebGraph[c1.front().m_source].location.x() << ", "
                      << reebGraph[c1.front().m_source].location.y() << "---" << std::endl;
            std::cout << "index : " << c1.back().m_target << " coords :"
                      << reebGraph[c1.back().m_target].location.x() << ", "
                      << reebGraph[c1.back().m_target].location.y() << "---" << std::endl;

            std::cout << "    chain2 :" << std::endl;
            std::cout << "index : " << c2.front().m_source << " coords : "
                      << reebGraph[c2.front().m_source].location.x() << ", "
                      << reebGraph[c2.front().m_source].location.y() << "---" << std::endl;
            std::cout << "index : " << c2.back().m_target << " coords : "
                      << reebGraph[c2.back().m_target].location.x() << ", "
                      << reebGraph[c2.back().m_target].location.y() << "---" << std::endl;

            std::cout << "    chain3 :" << std::endl;
            std::cout << "index : " << c3.front().m_source << " coords : "
                      << reebGraph[c3.front().m_source].location.x() << ", "
                      << reebGraph[c3.front().m_source].location.y() << "---" << std::endl;
            std::cout << "index : " << c3.back().m_target << " coords :"
                      << reebGraph[c3.back().m_target].location.x() << ", "
                      << reebGraph[c3.back().m_target].location.y() << "---" << std::endl
                      << std::endl;

            std::cout << "Step 3 : find the correct stem" << std::endl << std::endl;

            auto chainElementStartingFrom2 = [](std::vector<edge_descriptor> chain,
                                                int startElement, size_t idx) {
                if (chain.front().m_source == startElement) {
                    if (chain.size() == idx)
                        return chain.back().m_target;

                    return chain[idx].m_source;
                } else {
                    if (chain.size() == idx)
                        return chain.front().m_source;

                    return chain[chain.size() - 1 - idx].m_target;
                }
            };

            std::vector<std::vector<edge_descriptor>> threeChains = {c1, c2, c3};
            int yStemChain = -1;
            for (auto &adjChainIdx : {0, 1, 2}) {
                //if I can be connected to either of them, I am the stem
                size_t nextVertexAlongAdjChain = chainElementStartingFrom2(
                        threeChains[adjChainIdx], yJunctionIdx,
                        1);
                //simple hack to avoid trusting vertices way too close
                if ((threeChains[adjChainIdx].size() > 1) &&
                    (reebGraph[nextVertexAlongAdjChain].location -
                     reebGraph[yJunctionIdx].location).squaredNorm() < 1e-2)
                    nextVertexAlongAdjChain = chainElementStartingFrom2(threeChains[adjChainIdx],
                                                                        yJunctionIdx, 2);

                bool continuingBoth = true;
                for (auto &adjChainIdx2 : {0, 1, 2}) {
                    if (threeChains[adjChainIdx2] == threeChains[adjChainIdx])
                        continue;

                    size_t nextVertexAlongAdjChain2 = chainElementStartingFrom2(
                            threeChains[adjChainIdx2], yJunctionIdx, 1);
                    if ((threeChains[adjChainIdx2].size() > 1) &&
                        ((reebGraph[nextVertexAlongAdjChain2].location -
                          reebGraph[yJunctionIdx].location).squaredNorm() < 1e-2))
                        nextVertexAlongAdjChain2 = chainElementStartingFrom2(
                                threeChains[adjChainIdx2], yJunctionIdx, 2);

                    if (!yJunctionTest({nextVertexAlongAdjChain, (vertex_descriptor)yJunctionIdx},
                                       nextVertexAlongAdjChain2, reebGraph)) {
                        continuingBoth = false;
                        break;
                    }
                }
                if (continuingBoth) {
                    yStemChain = adjChainIdx;
                }
            }

            std::cout << "End of step 3 : choosing the correct stem. Correct stem : " << yStemChain
                      << std::endl << std::endl;

            if (yStemChain == -1) {
                //std::cout << "Can't find stem around vertex " << yStemChain << std::endl;
                continue;
            }

            std::cout << "Step 4 : Connecting the 3 chains into 2 chains sharing the stem"
                      << std::endl;

            //modifiedNewEdge = true;
            std::cout << "Found stem: chain : " << yStemChain << " out of 3. With the following next vertex : "
                      << chainElementStartingFrom2(threeChains[yStemChain], yJunctionIdx, 1)
                      << ". And y junction vertex : "<<yJunctionIdx<< std::endl;

            //now connect the rest to this chain
            std::map<size_t, std::vector<edge_descriptor>> newChainsMap;
            for (size_t adjChainIdx : {0, 1, 2}) {
                if (adjChainIdx == yStemChain)
                    continue;
                newChainsMap[adjChainIdx] = mergeChains(threeChains[adjChainIdx],
                                                        threeChains[yStemChain],
                                                        reebGraph);
            }

            if (yStemChain == 0) {
                chains[cIdx] = newChainsMap[1];
            } else if (yStemChain == 1) {
                chains[cIdx] = newChainsMap[0];
            } else if (yStemChain == 2) {
                chains[cIdx] = newChainsMap[0];
            }
            std::cout << "Now chains are " << std::endl;
            for (int i = 0; i<chains.size(); i++)
            {
                std::cout << "chain : " << i ;
                printChain(chains[i], reebGraph, true);
            }
        }
    }
//                    addCoverage(reebGraph, chainToVerticesSeq(newChain), TotalPixelCoverage, origMask);
//                    chains.push_back(newChain);
};

/**
 * This function removes adjacent chains that do not cover much. If such chain was found, we check that the other end is a labeled point.
 * If so, we remove the chain from the list, mark junction vertex as active,
 * and if (move_to_labeled_location) we move junction point to the labeled point location.
 * @param reebGraph
 * @param selectedJunctionPointIdx
 * @param selectedJunctionAdjChainsIdxs
 * @param selectedComponentId
 * @param activePtsPerComponent
 * @param chains
 * @param move_to_labeled_location
 * @param coverage_threshold
 * @return
 */
bool removeRedundantChains(
        G& reebGraph,
        const vertex_descriptor selectedJunctionPointIdx,
        const std::set<size_t>& selectedJunctionAdjChainsIdxs,
        const int selectedComponentId,
        std::vector<std::set<std::pair<int, int>>>& activePtsPerComponent,
        std::vector<std::vector<edge_descriptor>> &chains,
        bool move_to_labeled_location=true,
        int coverage_threshold=MINIMUM_COVERAGE_TO_BE_TAKEN_SERIOUSLY
        ) {

    bool fixedByMoving = false;
    std::vector<size_t> chainsToRemove;
    std::vector<size_t> vectorSelectedJunctionAdjChainsIdxs (selectedJunctionAdjChainsIdxs.begin(), selectedJunctionAdjChainsIdxs.end());
    // sorting the vector by coverage increasing
//                            std::sort(vectorSelectedJunctionAdjChainsIdxs.begin(), vectorSelectedJunctionAdjChainsIdxs.end(),
//                                 [&selectedJunctionPointIdx, &estimateCoverageOfAChainStartingFrom](const size_t & a, const size_t & b) -> bool
//                                 {
//                                     return estimateCoverageOfAChainStartingFrom(a, selectedJunctionPointIdx) < estimateCoverageOfAChainStartingFrom(b, selectedJunctionPointIdx);
//                                 });
    for (size_t adjChainIdx : vectorSelectedJunctionAdjChainsIdxs){
        // estimate coverage by this chain
        double coverageOfThisBranch = estimateCoverageOfAChainStartingFrom(reebGraph, chains, adjChainIdx, selectedJunctionPointIdx);
        vertex_descriptor otherEndpoint = getOtherEndOfChain(chains, adjChainIdx, selectedJunctionPointIdx);
        const int valenceOfOtherEndpoint = getVertexValence(otherEndpoint, chains);
        std::cout << "Coverage of a chain " << adjChainIdx << " is " << coverageOfThisBranch
                  << "; other endpoint is " << otherEndpoint << " with valence " << valenceOfOtherEndpoint << std::endl;
        // < threshold?
        if (coverageOfThisBranch < coverage_threshold) {
            if (fixedByMoving) {
                // if the chain is already solved: (1) remove the chain
                int labeledPointToInsert = getTheLabeledPointConnectedToMe(activePtsPerComponent, otherEndpoint, selectedComponentId);
                if (labeledPointToInsert>=0 && valenceOfOtherEndpoint==1) {
                    std::cout << "valence of endpoint " << otherEndpoint << " = " << getVertexValence(otherEndpoint, chains) << "\n";
                    chainsToRemove.push_back(adjChainIdx);
                }
            } else {
                // if the chain was not already solved: (1) move labeled point here (2) remove the chain
                int labeledPointToInsert = getTheLabeledPointConnectedToMe(activePtsPerComponent, otherEndpoint, selectedComponentId);
                if (labeledPointToInsert>=0) {
                    std::cout << "Marked this vertex (" << selectedJunctionPointIdx
                              << ") as active, as it is close to vertex " << otherEndpoint << " -- labeled point "
                              << labeledPointToInsert << "\n";
                    activePtsPerComponent[selectedComponentId].insert({ labeledPointToInsert, selectedJunctionPointIdx });
//                                        reebGraph[selectedJunctionPointIdx].location = reebGraph[otherEndpoint].location;
                    if (valenceOfOtherEndpoint==1) {
                        std::cout << "I should remove chain " << adjChainIdx << std::endl;
                        chainsToRemove.push_back(adjChainIdx);
                    } else {
                        std::cout << "I should remove chain " << adjChainIdx << std::endl;
                        chainsToRemove.push_back(adjChainIdx);
                        if (move_to_labeled_location)
                            reebGraph[selectedJunctionPointIdx].location = reebGraph[otherEndpoint].location;
                    }
                    fixedByMoving = true;
                }
            }
        } else {
            // since we sorted the vector by coverage increasing
//                                    break;
        }
    }
    if (fixedByMoving) {
        std::cout << "This junction is solved by moving!";
        std::sort(chainsToRemove.begin(), chainsToRemove.end(), std::greater<>());

        // if it is an unlabeled point - merge ramaining chains
        int labeledPointConnected = getTheLabeledPointConnectedToMe(activePtsPerComponent, selectedJunctionPointIdx, selectedComponentId);
        if ((labeledPointConnected == -1) && (chainsToRemove.size() == 1)) {
            //TODO: implement merging when you deleted just one adjacent chain
            std::cout << "Selected point is not labeled " << selectedJunctionPointIdx
                      << ", I should merge chains, IMPLEMENT ME!" << std::endl;
        }

        // erase unused chains and polys
        for (auto chainToRemoveIdx : chainsToRemove) {
            std::cout << "Removing chain " << chainToRemoveIdx << "\n";
            chains.erase(chains.begin() + chainToRemoveIdx);
        }
        std::cout << "New chains are " << std::endl;
        for (int i = 0; i < chains.size(); i++) {
            std::cout << "chain : " << i;
            printChain(chains[i], reebGraph, true);
        }
    }
    return fixedByMoving;
}

bool checkYJunctionBetweenChains(
        const G& reebGraph,
        const std::vector<edge_descriptor>& chain1,
        const std::vector<edge_descriptor>& chain2,
        const vertex_descriptor selectedJunctionPointIdx
){
    size_t nextVertexAlongAdjChain = chainElementStartingFrom(chain1, selectedJunctionPointIdx,1);
    if ((chain1.size()>1) && (reebGraph[nextVertexAlongAdjChain].location-reebGraph[selectedJunctionPointIdx].location).squaredNorm()<1e-2) // WN : If the vertex is too close to the start of the chain, take the second one insead. Shouldn't we check again until we are ok ? Are there some edge cases here ? Do we care ?
        nextVertexAlongAdjChain = chainElementStartingFrom(chain1, selectedJunctionPointIdx, 2);

    size_t nextVertexAlongAdjChain2 = chainElementStartingFrom(chain2, selectedJunctionPointIdx, 1);
    if ((chain2.size() > 1) && ((reebGraph[nextVertexAlongAdjChain2].location - reebGraph[selectedJunctionPointIdx].location).squaredNorm() < 1e-2))
        nextVertexAlongAdjChain2 = chainElementStartingFrom(chain2, selectedJunctionPointIdx, 2);

    return yJunctionTest({ nextVertexAlongAdjChain,selectedJunctionPointIdx }, nextVertexAlongAdjChain2, reebGraph);
}

std::pair<bool, int> findStem(
        const G& reebGraph,
        const std::vector<std::set<std::pair<int, int>>>& activePtsPerComponent,
        const int selectedComponentId,
        const vertex_descriptor selectedJunctionPointIdx,
        const std::set<size_t>& modifiedAdjChainsIdx,
        const std::vector<std::vector<edge_descriptor>> &chains
        ){
    int yStemChain = -1;
    for (size_t stemCandidateChainId : modifiedAdjChainsIdx){
        std::vector<edge_descriptor> stemCandidateChain = chains[stemCandidateChainId];
        bool canBeConnectedToBoth = true;
        for (size_t otherAdjacentChainId : modifiedAdjChainsIdx) {
            if (otherAdjacentChainId == stemCandidateChainId)
                continue;
            std::vector<edge_descriptor> otherAdjacentChain = chains[otherAdjacentChainId];
            canBeConnectedToBoth = canBeConnectedToBoth && checkYJunctionBetweenChains(reebGraph, stemCandidateChain, otherAdjacentChain, selectedJunctionPointIdx);
        }
        if (canBeConnectedToBoth) {
            std::cout << " OK WE FOUND A STEM CHAIN, IT IS : " << stemCandidateChainId << std::endl;
            yStemChain = stemCandidateChainId; // WN : Could this overwrite itself, or it cant happen more than once ?
            return std::make_pair(true, yStemChain);
        }
    }
    return std::make_pair(false, yStemChain);
}

/**
 * This function solves for all non-2 valence non-active junction points,
 * glue chains to form y-junction (with overlapping)
 * @param reebGraph - reebGraph, is changed only if redundant lines are pruned.
 * @param pts
 * @param selectedComponentId
 * @param chains
 * @param activePtsPerComponent
 */
void solveJunctions(
        G& reebGraph,
        const std::vector<cv::Point2d>& pts,
        const int selectedComponentId,
        std::vector<std::vector<edge_descriptor>>& chains,
        std::vector<std::set<std::pair<int, int>>>& activePtsPerComponent
        ){
    std::map<size_t, std::set<size_t>> unexpectedEndChains_adjacentChains; // map from vertex id to connected chains (only non-active points here)
    std::vector<size_t> stemCompensation(boost::num_vertices(reebGraph), 0); // vector of stem compensation indices

    do {
        unexpectedEndChains_adjacentChains.clear();

        unexpectedEndChains_adjacentChains = findUnexpectedEndChainsVertices(reebGraph, pts, chains,
                                                                             activePtsPerComponent,
                                                                             selectedComponentId);

        if (!unexpectedEndChains_adjacentChains.empty()) {
            auto firstelem = *unexpectedEndChains_adjacentChains.begin();
            vertex_descriptor selectedJunctionPointIdx = firstelem.first; // id of a junction point we process
            std::set<size_t> selectedJunctionAdjChainsIdxs = firstelem.second; // set (ids of adjacent chains)
            int selectedJunctionPointValence = selectedJunctionAdjChainsIdxs.size(); // valence of the junction point

            // print some info
            std::cout << "I should remove vertex " << selectedJunctionPointIdx << " valence: "
                      << selectedJunctionPointValence << " (with stem compensation = "
                      << stemCompensation[selectedJunctionPointIdx] << ")" << std::endl;
//                // If you don't want to solve for junction with > 3 - use this
//            if (selectedJunctionPointValence>3){
//                std::cout << "I'm not touching vertex " << selectedJunctionPointIdx << " at component "
//                          << selectedComponentId << " (# of adjacent chains: "
//                          << selectedJunctionPointValence << ")" << std::endl;
//
//                activePtsPerComponent[selectedComponentId].insert({-1,
//                                                                   selectedJunctionPointIdx});
//                continue;
//
//            }

            for (auto e : selectedJunctionAdjChainsIdxs) {
                std::cout << "adjacent chain : " << e;
                printChain(chains[e], reebGraph, true);
            }

            // check if we can fix this junction point if it is close to the labeled point
            bool fixedByMoving = removeRedundantChains(reebGraph, selectedJunctionPointIdx,
                                                       selectedJunctionAdjChainsIdxs,
                                                       selectedComponentId, activePtsPerComponent,
                                                       chains);

            // if not, solve the junction ()
            if (selectedJunctionPointValence - stemCompensation[selectedJunctionPointIdx] == 3 &&
                stemCompensation[selectedJunctionPointIdx] <= 1 &&
                !fixedByMoving) // We care only if valence == 3. ADDING A TEST FOR VALENCE > 3 CAUSED BY DOUBLE EDGES WIth STEMS
            {
                // print info
                std::cout << "WORKING WITH VERTEX :" << selectedJunctionPointIdx << ". Location : "
                          << reebGraph[selectedJunctionPointIdx].location.x() << ", "
                          << reebGraph[selectedJunctionPointIdx].location.y() << std::endl;
                std::cout << "chains in work: ";
                for (auto cc : selectedJunctionAdjChainsIdxs)
                    std::cout << cc << " ";
                std::cout << std::endl;

                size_t yStemChain = -1;

                // The idea here is to remove duplicated chain (happens when you simplify valence 3 junction) and another valence 3 happens near
                // Playing around with the adjacent chains idx so that we dont search for stem with doubled chains (happens when stemCompensation[selectedJunctionPointIdx] is 1.
                std::set originalAdjChainsIdx = selectedJunctionAdjChainsIdxs;
                std::set modifiedAdjChainsIdx = selectedJunctionAdjChainsIdxs;
                bool foundChainToIgnore = false;
                std::vector<edge_descriptor> c1, c2;
                size_t previousStem1; // index of the first previous stem
                size_t previousStem2; // index of the second previous stem
                for (size_t i : originalAdjChainsIdx) {
                    c1 = chains[i];

                    for (auto j : originalAdjChainsIdx) {
                        if (i == j) continue;
                        c2 = chains[j];
                        for (auto e1 : c1) {
                            for (auto e2 : c2) {
                                if (!foundChainToIgnore && e1.m_source == e2.m_source &&
                                    e1.m_source != selectedJunctionPointIdx) {
                                    auto it_chainToIgnore = modifiedAdjChainsIdx.find(i);
                                    modifiedAdjChainsIdx.erase(it_chainToIgnore);
                                    foundChainToIgnore = true;
                                    previousStem1 = i;
                                    previousStem2 = j;
                                    std::cout << "Found chain to delete: " << i << " overlaps with "
                                              << j << ", remove the former one for some time"
                                              << std::endl;
                                    break;
                                }
                            }
                        }
                    }
                }

                std::pair<bool, int> stemSearchResult = findStem(reebGraph, activePtsPerComponent, selectedComponentId,
                                      selectedJunctionPointIdx, modifiedAdjChainsIdx, chains);

                if (!stemSearchResult.first) {
                    std::cout << "Can't find stem around vertex " << yStemChain << std::endl;

                } else {
                    yStemChain = stemSearchResult.second;
                    std::cout << "Found stem: " << yStemChain
                              << " chain with the following next vertex "
                              << chainElementStartingFrom(chains, yStemChain, selectedJunctionPointIdx,
                                                          1) << std::endl;
                    printChain(chains[yStemChain], reebGraph, true);
                }

                std::map<size_t, std::vector<edge_descriptor>> newChainAndPolyline;
                //now connect the rest to the stem chain
                if (
                        (stemCompensation[selectedJunctionPointIdx] == 0) ||
                        // normal case, junction of valence 3
                        (stemCompensation[selectedJunctionPointIdx] == 1 &&
                         yStemChain != previousStem1 && yStemChain !=
                                                        previousStem2) // valence 4 and new stem is not the old stem (from some previous y-junction merge)
                        ) {
                    // merge every other chain with the stemChain
                    for (size_t adjChainIdx : originalAdjChainsIdx) {
                        if (adjChainIdx == yStemChain)
                            continue;
                        newChainAndPolyline[adjChainIdx] = mergeChains(chains[adjChainIdx],
                                                                       chains[yStemChain], reebGraph);
                        // do the merging
                        chains[adjChainIdx] = mergeChains(chains[adjChainIdx], chains[yStemChain],
                                                          reebGraph);

                    }
                    // now remove the stem chain
                    chains.erase(chains.begin() + yStemChain);
                } else if (stemCompensation[selectedJunctionPointIdx] ==
                           1) // valence 4 and new stem is one of the previous y-junctions
                {
                    if (previousStem2 == yStemChain) {
                        bool mergedOne = false;
                        for (size_t adjChainIdx : modifiedAdjChainsIdx) // WN : Iterate over the chains adjacent to vertex selectedJunctionPointIdx
                        { // WN : Merge each chain that is not the yStemChain to the yStemchain.
                            if (adjChainIdx == previousStem1 || adjChainIdx == previousStem2)
                                continue;

                            int mergeStem = mergedOne ? previousStem2 : previousStem1;
                            newChainAndPolyline[adjChainIdx] = mergeChains(chains[adjChainIdx],
                                                                           chains[mergeStem],
                                                                           reebGraph);
                            mergedOne = true;
                            //printChain(newChainAndPolyline[adjChainIdx].first, reebGraph, false);
                        }
                        for (size_t adjChainIdx : modifiedAdjChainsIdx) // WN : Iterate over the chains adjacent to vertex selectedJunctionPointIdx
                        {
                            if (adjChainIdx == previousStem1 || adjChainIdx == previousStem2)
                                continue;

                            chains[adjChainIdx] = newChainAndPolyline[adjChainIdx];
                        }
                        chains.erase(chains.begin() + previousStem1);

                        previousStem2 =
                                previousStem1 < previousStem2 ? previousStem2 - 1 : previousStem2;

                        chains.erase(chains.begin() + previousStem2);
                        stemCompensation[selectedJunctionPointIdx]--;
                    } else
                        std::cout << "Stem chain is not one of the previous stem chains" << std::endl;
                }


                // I wanna mark the end of the stem that isnt the y junction
                std::vector<size_t> endPts;
                for (size_t adjChainIdx : modifiedAdjChainsIdx) {
                    if (adjChainIdx == yStemChain)
                        continue;
                    endPts.push_back(newChainAndPolyline[adjChainIdx].front().m_source);
                    endPts.push_back(newChainAndPolyline[adjChainIdx].back().m_target);
                }

                std::map<size_t, int> countMap;
                // Iterate over the vector and store the frequency of each element in map
                for (auto &elem : endPts) {
                    auto result = countMap.insert(std::pair<size_t, int>(elem, 1));
                    if (result.second == false)
                        result.first->second++;
                }

                for (auto elem : countMap) {
                    if (elem.second >= 2)
                        stemCompensation[elem.first] += elem.second - 1;
                }

                std::cout << "Removed vertex " << selectedJunctionPointIdx << std::endl;
                std::cout << "New chains are " << std::endl;
                for (int i = 0; i < chains.size(); i++) {
                    std::cout << "chain : " << i;
                    printChain(chains[i], reebGraph, true);
                }

            } else if (selectedJunctionPointValence >=
                       5 /*&& stemCompensation[selectedJunctionPointIdx] >= 0*/ ) {
                stemCompensation[selectedJunctionPointIdx] = -5; // ?????
                std::cout << " ** VALENCE 5 or more ** WORKING WITH VERTEX :"
                          << selectedJunctionPointIdx << ". Location : "
                          << reebGraph[selectedJunctionPointIdx].location.x() << ", "
                          << reebGraph[selectedJunctionPointIdx].location.y()
                          << std::endl;
                for (int cIdx : selectedJunctionAdjChainsIdxs) {
                    std::cout << "Chain in val 5 : " << cIdx;
                    printChain(chains[cIdx], reebGraph, true);
                    std::cout << "chain length : " << chainLength(chains[cIdx], reebGraph) << std::endl
                              << std::endl;
                }

                std::set<std::pair<size_t, size_t>> pairChains;

                //take 2 branches at once. For now just find pairs and print their id
                for (int cIdx1 : selectedJunctionAdjChainsIdxs) {
                    size_t nextVertexAlongAdjChain = chainElementStartingFrom(chains, cIdx1,
                                                                              selectedJunctionPointIdx,
                                                                              1); // WN : Next vertex along the chain
                    if ((chains[cIdx1].size() > 1) && (reebGraph[nextVertexAlongAdjChain].location -
                                                       reebGraph[selectedJunctionPointIdx].location).squaredNorm() <
                                                      1e-2) // WN : If the vertex is too close to the start of the chain, take the second one insead. Shouldn't we check again until we are ok ? Are there some edge cases here ? Do we care ?
                        nextVertexAlongAdjChain = chainElementStartingFrom(chains, cIdx1,
                                                                           selectedJunctionPointIdx, 2);
                    for (int cIdx2 : selectedJunctionAdjChainsIdxs) {
                        if (cIdx1 == cIdx2) continue;
                        size_t nextVertexAlongAdjChain2 = chainElementStartingFrom(chains, cIdx2,
                                                                                   selectedJunctionPointIdx,
                                                                                   1);
                        if ((chains[cIdx2].size() > 1) &&
                            ((reebGraph[nextVertexAlongAdjChain2].location -
                              reebGraph[selectedJunctionPointIdx].location).squaredNorm() < 1e-2))
                            nextVertexAlongAdjChain2 = chainElementStartingFrom(chains, cIdx2,
                                                                                selectedJunctionPointIdx,
                                                                                2);

                        if (yJunctionTest({nextVertexAlongAdjChain, selectedJunctionPointIdx},
                                          nextVertexAlongAdjChain2,
                                          reebGraph)) // WN : selectedJunctionPointIdx is the annotation point. nextVertexAlongAdjChain is it's neighbor on one of the adjacent chain. nextVertexAlongAdjChain2 is another neighbor on another chain Todo read yJunctionTest()
                        { // WN : If this returns false, stop there. The current value of yStemChain will determine if and which chain will become a stem for its neighbors.
                            std::cout << "should be guuuud between chain " << cIdx1 << " + " << cIdx2
                                      << std::endl;
                            auto cMin = cIdx1 < cIdx2 ? cIdx1 : cIdx2;
                            auto cMax = cIdx1 < cIdx2 ? cIdx2 : cIdx1;
                            pairChains.insert(std::pair(cMin, cMax));
                        }
                    }
                }

                std::vector<std::pair<std::pair<size_t, size_t>, std::vector<edge_descriptor>>> tentativeMergedChains;
                for (auto twoChains : pairChains) {
                    std::cout << twoChains.first << " " << twoChains.second << std::endl;
                    //MyPolyline emptyPoly1, emptyPoly2;
                    auto newChain = mergeChains(chains[twoChains.first], chains[twoChains.second],
                                                reebGraph);
                    tentativeMergedChains.push_back(std::pair(twoChains, newChain));
                }

                std::vector<std::pair<std::pair<size_t, size_t>, std::vector<edge_descriptor>>> in;
                std::vector<std::vector<edge_descriptor>> finalChoice;
                if (checkChainOfChains(in, tentativeMergedChains, selectedJunctionAdjChainsIdxs,
                                       finalChoice, reebGraph)) {
                    std::cout << "wow found smth !!!" << std::endl;
                    for (auto e : finalChoice)
                        printChain(e, reebGraph, true);
                }

                for (auto cIdx : selectedJunctionAdjChainsIdxs) {
                    if (finalChoice.size() != 0) {
                        chains[cIdx] = finalChoice[0];
                        finalChoice.erase(finalChoice.begin());
                    } else {
                        chains.erase(chains.begin() + cIdx);
                    }
                }

            } else // WN : when valence is not 3.
            {
                std::cout << "I'm not touching vertex " << selectedJunctionPointIdx << " at component "
                          << selectedComponentId << " (# of adjacent chains: "
                          << selectedJunctionPointValence << ")" << std::endl;
            }

            activePtsPerComponent[selectedComponentId].insert({-1,
                                                               selectedJunctionPointIdx}); //mark that we're not touching this vertex now. WN : Reminder, selectedJunctionPointIdx is the 'annotation point' (aka an end of the chain close enough to a labelled point).
        }
    } while (!unexpectedEndChains_adjacentChains.empty());
}

/**
 * This funtion checks every unexpected junction and glues together all chains that pass yJunctionTest
 * @param reebGraph
 * @param pts
 * @param selectedComponentId
 * @param chains
 * @param activePtsPerComponent
 */
void glueAllJunctions(
        G& reebGraph,
        const std::vector<cv::Point2d>& pts,
        const int selectedComponentId,
        std::vector<std::vector<edge_descriptor>>& chains,
        std::vector<std::set<std::pair<int, int>>>& activePtsPerComponent
){
    std::map<size_t, std::set<size_t>> unexpectedEndChains_adjacentChains; // map from vertex id to connected chains (only non-active points here)

    do {
        unexpectedEndChains_adjacentChains.clear();

        unexpectedEndChains_adjacentChains = findUnexpectedEndChainsVertices(reebGraph, pts, chains,
                                                                             activePtsPerComponent,
                                                                             selectedComponentId);

        if (!unexpectedEndChains_adjacentChains.empty()) {
            auto firstelem = *unexpectedEndChains_adjacentChains.begin();
            vertex_descriptor selectedJunctionPointIdx = firstelem.first; // id of a junction point we process
            std::set<size_t> selectedJunctionAdjChainsIdxs = firstelem.second; // set (ids of adjacent chains)
            std::vector<size_t> selectedJunctionAdjChainsIdxsVector (selectedJunctionAdjChainsIdxs.begin(), selectedJunctionAdjChainsIdxs.end());
            // sort indices in descending order (useful for removing)
            std::sort(selectedJunctionAdjChainsIdxsVector.begin(), selectedJunctionAdjChainsIdxsVector.end(), std::greater<>());
            int selectedJunctionPointValence = selectedJunctionAdjChainsIdxs.size(); // valence of the junction point

            // print some info
            std::cout << "I found unexpected vertex " << selectedJunctionPointIdx << " valence: "
                      << selectedJunctionPointValence << std::endl;

            for (auto e : selectedJunctionAdjChainsIdxs) {
                std::cout << "adjacent chain : " << e;
                printChain(chains[e], reebGraph, true);
            }

            if (selectedJunctionPointValence==1){
                std::cout << "Valence 1, let it live.\n";
                activePtsPerComponent[selectedComponentId].insert({-1, selectedJunctionPointIdx});
                continue;
            }

            std::set<std::pair<size_t, size_t>> pairsToGlue;
            for (size_t chainInAdjacent_ListId=0; chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++chainInAdjacent_ListId){
                size_t basechainId = selectedJunctionAdjChainsIdxsVector[chainInAdjacent_ListId];
                std::cout << "working with chain: " << basechainId << "\n";
                Chain basechain = chains[basechainId];
                for (size_t other_chainInAdjacent_ListId=chainInAdjacent_ListId+1; other_chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++other_chainInAdjacent_ListId){
                    size_t chainToGlueWithId = selectedJunctionAdjChainsIdxsVector[other_chainInAdjacent_ListId];
                    Chain chainToGlueWith = chains[chainToGlueWithId];
                    bool yjunctonIsCorrect = checkYJunctionBetweenChains(reebGraph, basechain, chainToGlueWith, selectedJunctionPointIdx);
                    if (yjunctonIsCorrect) {
                        std::cout<<"Yes, chains " << basechainId << " and " << chainToGlueWithId << " can be connected\n";
                        pairsToGlue.insert(std::make_pair(basechainId, chainToGlueWithId));
                    } else {
                        std::cout<<"No, I can't connect chains " << basechainId << " and " << chainToGlueWithId << "\n";
                    }
                }

            }
            assert(!pairsToGlue.empty());

            // now we glue all plausible pairs
            std::vector<Chain> newChains;
            for (auto pairToGlue : pairsToGlue){
                newChains.push_back(mergeChains(chains[pairToGlue.first], chains[pairToGlue.second], reebGraph));
            }
            // remove old chains
            for (size_t chainToRemove : selectedJunctionAdjChainsIdxsVector) {
                chains.erase(chains.begin() + chainToRemove);
            }
            // add new chains
            for (auto newchain : newChains) {
                chains.push_back(newchain);
            }

            std::cout << "Now chains are " << std::endl;
            for (int i = 0; i < chains.size(); i++) {
                std::cout << "chain : " << i;
                printChain(chains[i], reebGraph, true);
            }

        }

    } while (!unexpectedEndChains_adjacentChains.empty());
}

std::vector<Chain> advancedGlueAllJunction (
        G& reebGraph,
        const std::vector<cv::Point2d>& pts,
        const int selectedComponentId,
        std::vector<std::vector<edge_descriptor>>& chains,
        std::vector<std::set<std::pair<int, int>>>& activePtsPerComponent
    ){
    std::map<size_t, std::set<size_t>> unexpectedEndChains_adjacentChains; // map from vertex id to connected chains (only non-active points here)

    do {
        unexpectedEndChains_adjacentChains.clear();

        unexpectedEndChains_adjacentChains = findUnexpectedEndChainsVertices(reebGraph, pts, chains,
                                                                             activePtsPerComponent,
                                                                             selectedComponentId);

        if (!unexpectedEndChains_adjacentChains.empty()) {
            auto firstelem = *unexpectedEndChains_adjacentChains.begin();
            vertex_descriptor selectedJunctionPointIdx = firstelem.first; // id of a junction point we process
            std::set<size_t> selectedJunctionAdjChainsIdxs = firstelem.second; // set (ids of adjacent chains)
            std::vector<size_t> selectedJunctionAdjChainsIdxsVector (selectedJunctionAdjChainsIdxs.begin(), selectedJunctionAdjChainsIdxs.end());
            // sort indices in descending order (useful for removing)
            std::sort(selectedJunctionAdjChainsIdxsVector.begin(), selectedJunctionAdjChainsIdxsVector.end(), std::greater<>());
            int selectedJunctionPointValence = selectedJunctionAdjChainsIdxs.size(); // valence of the junction point

            // print some info
            std::cout << "I found unexpected vertex " << selectedJunctionPointIdx << " valence: "
                      << selectedJunctionPointValence << std::endl;

            for (auto e : selectedJunctionAdjChainsIdxs) {
                std::cout << "adjacent chain : " << e;
                printChain(chains[e], reebGraph, true);
            }

            if (selectedJunctionPointValence==1){
                std::cout << "Valence 1, let it live.\n";
                activePtsPerComponent[selectedComponentId].insert({-1, selectedJunctionPointIdx});
                continue;
            }

            std::vector<Chain> chainInAdjacent_List;
            for (size_t chainInAdjacent_ListId=0; chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++chainInAdjacent_ListId){
                Chain chainToPush = chains[selectedJunctionAdjChainsIdxsVector[chainInAdjacent_ListId]];
                chainToPush = alignChainToStartWith(reebGraph, selectedJunctionPointIdx, chainToPush);
                chainInAdjacent_List.push_back(chainToPush);
                for (size_t otherChain_ListId=0; otherChain_ListId<chainInAdjacent_ListId; ++otherChain_ListId){
                    Chain chainPushed = chainInAdjacent_List[chainInAdjacent_ListId];
                    Chain otherChain = chainInAdjacent_List[otherChain_ListId];
                    if (chainsOverlapFromStart(chainPushed, otherChain)){
                        std::cout << "Jeez, we overlap! " << selectedJunctionAdjChainsIdxsVector[chainInAdjacent_ListId] << " and " << selectedJunctionAdjChainsIdxsVector[otherChain_ListId] << "\n";
                        std::pair<Chain, Chain> chain1decomposed, chain2decomposed; // pair (overlap, residual)
                        std::tie(chain1decomposed, chain2decomposed) = decomposeOverlappingChains(chainPushed, otherChain);
                        double overlapLength = chainLength(chain1decomposed.first, reebGraph);
                        std::cout << "Overlap length: " << overlapLength << "\n";
                        if (overlapLength<10){
                            std::cout << "Overlap length lower than threshold, cropping " << overlapLength << "\n";
                            chainInAdjacent_List[chainInAdjacent_ListId] = chain1decomposed.second;
                            chainInAdjacent_List[otherChain_ListId] = chain2decomposed.second;
                        }

                    }
                }
            }

            for (size_t chainInAdjacent_ListId=0; chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++chainInAdjacent_ListId){
                std::cout << "(flipped, nonoverlap) adjacent chain : " << selectedJunctionAdjChainsIdxsVector[chainInAdjacent_ListId];
                printChain(chainInAdjacent_List[chainInAdjacent_ListId], reebGraph, true);
            }

            std::set<std::pair<size_t, size_t>> pairsToGlue;
            for (size_t chainInAdjacent_ListId=0; chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++chainInAdjacent_ListId){
                size_t basechainId = selectedJunctionAdjChainsIdxsVector[chainInAdjacent_ListId];
                std::cout << "working with chain: " << basechainId << "\n";
                Chain basechain = chains[basechainId];
                Chain basechain_nonoverlap = chainInAdjacent_List[chainInAdjacent_ListId];
                for (size_t other_chainInAdjacent_ListId=chainInAdjacent_ListId+1; other_chainInAdjacent_ListId<selectedJunctionAdjChainsIdxsVector.size(); ++other_chainInAdjacent_ListId){
                    size_t chainToGlueWithId = selectedJunctionAdjChainsIdxsVector[other_chainInAdjacent_ListId];
                    Chain chainToGlueWith = chains[chainToGlueWithId];
                    Chain chainToGlueWith_nonoverlap = chainInAdjacent_List[other_chainInAdjacent_ListId];

                    size_t nextVertexAlongBaseChain = basechain_nonoverlap.front().m_target;
                    size_t nextVertexAlongChainToGlue = chainToGlueWith_nonoverlap.front().m_target;
                    std::cout << "Taking vertex " << nextVertexAlongBaseChain << " and " << nextVertexAlongChainToGlue << "\n";
//                    bool yjunctonIsCorrect = checkYJunctionBetweenChains(reebGraph, basechain, chainToGlueWith, selectedJunctionPointIdx);
                    bool yjunctonIsCorrect = yJunctionTest({nextVertexAlongBaseChain, selectedJunctionPointIdx}, nextVertexAlongChainToGlue, reebGraph);
                    if (yjunctonIsCorrect) {
                        size_t curvesShared = countSharedCurves(nextVertexAlongBaseChain, nextVertexAlongChainToGlue, reebGraph);
                        std::cout<<"Yes, chains " << basechainId << " and " << chainToGlueWithId
                        << " can be connected ("<< curvesShared <<" curves shared)\n";
                        pairsToGlue.insert(std::make_pair(basechainId, chainToGlueWithId));
//                        // This condition is too harsh, discarding it
//                        if (curvesShared>0)
//                            pairsToGlue.insert(std::make_pair(basechainId, chainToGlueWithId));
//                        else
//                            std::cout << "But I will not, because no shared curves\n";
                    } else {
                        std::cout<<"No, I can't connect chains " << basechainId << " and " << chainToGlueWithId << "\n";
                    }
                }

            }
            assert(!pairsToGlue.empty());

            // now we glue all plausible pairs
            std::vector<Chain> newChains;
            for (auto pairToGlue : pairsToGlue){
                newChains.push_back(mergeChains(chains[pairToGlue.first], chains[pairToGlue.second], reebGraph));
            }
            // remove old chains
            for (size_t chainToRemove : selectedJunctionAdjChainsIdxsVector) {
                chains.erase(chains.begin() + chainToRemove);
            }
            // add new chains
            for (auto newchain : newChains) {
                chains.push_back(newchain);
            }

            std::cout << "Now chains are " << std::endl;
            for (int i = 0; i < chains.size(); i++) {
                std::cout << "chain : " << i;
                printChain(chains[i], reebGraph, true);
            }

        }

    } while (!unexpectedEndChains_adjacentChains.empty());
    return chains;
}

/**
 * Main function to: build Stein tree, add more connections maximizing coverage, solve junctions
 * @param reebGraph
 * @param lpts - coordinates of labeled points
 * @param bwImg
 * @param origMask - boolean mask of the image
 * @param covered - returns a map from vertices to boolean (covered / uncovered)
 * @param fff - framefield flow
 * @param allRoots
 * @param newlyAddedIntersections - return a vector of endpoints we added to increase coverage
 * @param finalchains - return a vector of (vector of vertices id) in chains
 * @param allowedChains - all allowed chains
 * @param keypointToValenceMap - labeled keypoint to its valence
 * @param allowedChainToLabeledPointPair - each allowed chain to its labeled points
 * @return
 */
std::vector<MyPolyline> optimizeTopology(
        G& reebGraph,
        const std::vector<cv::Point2d>& lpts,
        const cv::Mat& bwImg,
        const cv::Mat& origMask,
        std::map<vertex_descriptor,
        bool>& covered,
        const FrameFieldFlow& fff,
        const std::array<Eigen::MatrixXcd, 2>& allRoots,
        std::vector<cv::Point2d>& newlyAddedIntersections,
        std::vector<std::vector<vertex_descriptor>>& finalchains,
        std::vector<std::vector<vertex_descriptor>>& allowedChains,
        std::map<int, int>& keypointToValenceMap,
        std::vector<std::pair<int, int>>& allowedChainToLabeledPointPair
        )
{
    std::vector<cv::Point2d> pts = lpts; // labeled points coordinates
    std::vector<MyPolyline> polylineresult;
    std::vector<int> vertexToComponent; // vector from vertex to component id. vertexToComponent[i] is which component vertex i is a part of
    std::vector<std::vector<vertex_descriptor>> componentToVertices; // vector from component id to vector of vertices. componentToVertices[j] is which vertices are part of component j
    std::tie(vertexToComponent, componentToVertices) = computeComponents(reebGraph);
    int nComponents = componentToVertices.size(); // number of components

    finalchains.clear();
    allowedChains.clear();
    keypointToValenceMap.clear();
    allowedChainToLabeledPointPair.clear();

    std::vector<cv::Point2d> newFarthestPointsToAdd_coords; // contains coordinates of new active points we add to the graph
    std::vector<size_t> newFarthestPointsToAdd; // contains indicies of new active points we add to the  graph

    std::set<std::pair<int, int>> TotalPixelCoverage; // set of covered pixels coordinates

    //now, for each component let's find all active pts
    const double activePointDistanceThreshold = 0; // used to be 4
    std::vector<std::set<std::pair<int, int>>> activePtsPerComponent; // vector of active points. activePtsPerComponent[a] = {{b, c}, ...} where a = component, b = labelled pt, c = closest point to b of component a.
    activePtsPerComponent = computeActivePointsPerComponent(reebGraph, pts, origMask, nComponents, vertexToComponent, allRoots, false, false);

    std::queue<int> componentsToProcess; // queue of components id to process
//    std::set<size_t> allActiveVerts; // set of all active vertices

    //now add junctions at sharp corners where there is now label
    std::cout<<"Here the addJunctions is running (commented)"<<std::endl<<std::endl;
//    newlyAddedIntersections = addJunctions(reebGraph, allActiveVerts, vertexToComponent, pts, activePointDistanceThreshold); // WN : Fct returns a list of vertex coordinates, already added to g.

//    //if you did addJunctions - recompute the vertexToComponent and active pts
//    std::tie(vertexToComponent, componentToVertices) = computeComponents(reebGraph); // WN : The number of component stays the same afaik, just added vertices at intersections of 2 vertexToComponent.
//    activePtsPerComponent = computeActivePointsPerComponent(reebGraph, pts, origMask, nComponents, vertexToComponent, allRoots, true);
//    std::cout<<"Number of components : "<<nComponents<<std::endl<<std::endl;

    std::cout<<"----------------------\nPRINTING COMPONENTS\n----------------------\n";
    for (int component = 0; component < nComponents; ++component)
    {
        // discarding a component if it is too short (0 or 1 vertex in the component)
        if (componentToVertices[component].size()<2)
            continue;
        componentsToProcess.push(component);
        std::cout<< "COMPONENT: " << component <<".  Active vertices: ";
        for (const auto& a : activePtsPerComponent[component]) {
//            allActiveVerts.insert(a.second);
            std::cout << " ( labeled point "<<a.first<<", closest vertex "<<a.second<< ") - ";
        }
        std::cout << "\n All vertices : ";
        for (size_t v=0; v < boost::num_vertices(reebGraph); ++v) {
            if (vertexToComponent[v] == component)
                std::cout << v << ", ";
        }
        std::cout<< " - END COMPONENT" <<std::endl;
    }

    std::cout<<"----------------------\nPROCESSING COMPONENTS\n----------------------\n";
    do { // process each component in the queue
        int selectedComponentId = componentsToProcess.front();

        if ((activePtsPerComponent[selectedComponentId].size() > 1) && (activePtsPerComponent[selectedComponentId].size() < MAXIMUM_ACTIVE_POINTS)) {
            std::set<std::pair<int, int>> activePointsPairs = activePtsPerComponent[selectedComponentId];
            activePtsPerComponent[selectedComponentId] = filterActivePtsByDistance(reebGraph, activePointsPairs, pts, true);
        }

        // print active points in the component
        std::cout << "\n\n\n----------\n--> COMPONENT " << selectedComponentId << "\n End pts: ";
        for (auto v : activePtsPerComponent[selectedComponentId])
            std::cout << "\n" << v.first << "\t (graph vertex " << v.second << "),\t coords: " << reebGraph[v.second].location.transpose();
        std::cout << std::endl;

        std::set<size_t> activepoints; // set of active points of current component.
        for (const auto& it: activePtsPerComponent[selectedComponentId])
            activepoints.insert(it.second);

        //now for active vertices, if there are more than 1, let's find a Steiner tree
        if (activepoints.size() > 1)
        {

            std::set<vertex_descriptor> terminals;
            for (const auto &ac : activePtsPerComponent[selectedComponentId]) {
                terminals.insert(
                        ac.second); // WN : Reminder, ac.first is a labelled point, ac.second is the closest graph vertex to it.
            }

            std::map<vertex_descriptor, bool> vertexIsCoveredByRadius; // map, from vertex index -> true if covered (close to a edge from chain)
            bool coverSomeMore = true;
            if (terminals.size() < MAXIMUM_ACTIVE_POINTS) {

                // get chains from steiner tree on active vertices
                std::vector<std::vector<edge_descriptor>> chains = buildChainsFromSteinerTree(reebGraph, activepoints);

                // print chains
                for (int i = 0; i<chains.size(); i++)
                {
                    std::cout << "chain after steiner : " << i;
                    printChain(chains[i], reebGraph, true);
                }

                // print endpoints and update pixel coverage
                std::cout<<"START AND END VERTICES"<<std::endl;
                for (int cIdx = 0; cIdx<chains.size();cIdx++)
                {
                    std::cout<<chains[cIdx].front().m_source<<" --- "<<chains[cIdx].back().m_target<<std::endl;
                    addCoverage(reebGraph, chainToVerticesSeq(chains[cIdx]), TotalPixelCoverage, origMask);
                }



                //ok great, now we have Steiner trees. let's now extend them with shortest edges each covering at least a given number of pixels
                //I'm using an algorithm similar to https://link.springer.com/article/10.1007/s10707-019-00385-8

                for (const auto &c : chains) {
                    //std::cout << " chain " << std::endl;
                    markVerticesCoveredByRadius(reebGraph, c, componentToVertices, selectedComponentId, vertexIsCoveredByRadius);
                }
//                std::set<vertex_descriptor> terminals;
//                for (const auto &ac : activePtsPerComponent[selectedComponentId]) {
//                    terminals.insert(
//                            ac.second); // WN : Reminder, ac.first is a labelled point, ac.second is the closest graph vertex to it.
//                }
                for(auto t : terminals)
                    std::cout<<"terminal : "<<t<<std::endl;
                //now start adding shortest edges one by one, each covering at least N new vertices

                std::cout << "Looking for uncovered areas..." << std::endl; //TODO: optimize using the topograph
                std::vector<edge_descriptor> newEdge;

                std::cout << " - Trying to add new paths without adding new endpoints\n";
                do {
                    const int N = 4; // WN : Minimum amount of vertices to cover ?
                    newEdge = findShortestPath_AtLeastNNotCoveredVerticesMNewPixels(reebGraph, componentToVertices[selectedComponentId], terminals,
                                                                          vertexIsCoveredByRadius, N,
                                                                          fff,
                                                                          TotalPixelCoverage,
                                                                          MINIMUM_PIXEL_COVERAGE ,
                                                                          origMask);
                    markVerticesCoveredByRadius(reebGraph, newEdge, componentToVertices, selectedComponentId, vertexIsCoveredByRadius);

                    if (!newEdge.empty()) {
                        tryToAddNewChain(reebGraph, origMask, newEdge, chains, TotalPixelCoverage);
                    }
                    else {
                        std::cout << "new edge is empty!\n";
                    }
                } while (!newEdge.empty()); // WN : aka didnt find any new edges to add.

                std::cout << " - Adding new endpoints. Trying to add new paths.\n";
                do {
                    std::set<size_t> coveredVertices = getVerticesUsedInChains(chains);
                    auto shouldIAddNewPointForThisComponent = findFarthestUncoveredVertexFromTheSameComponent(reebGraph, coveredVertices, 7.0);
//                    auto shouldIAddNewPointForThisComponent = findFarthestPixelUncoveredVertexFromTheSameComponent(reebGraph, coveredVertices, TotalPixelCoverage, 2);
                    if (shouldIAddNewPointForThisComponent.first) {
                        size_t newpoint_index = shouldIAddNewPointForThisComponent.second;
                        std::cout << "Adding uncovered vertex " << newpoint_index << " as a terminal to component " << selectedComponentId << std::endl;
                        terminals.insert(newpoint_index);
                        newFarthestPointsToAdd.push_back(newpoint_index);
                        newFarthestPointsToAdd_coords.emplace_back( reebGraph[newpoint_index].location.y(), reebGraph[newpoint_index].location.x() );
                        activePtsPerComponent[selectedComponentId].insert(std::make_pair(-newFarthestPointsToAdd.size(), newpoint_index));
                    }
                    const int N = 10; // WN : Minimum amount of vertices to cover ?
                    //TODO: maybe do not use findShortestPath_AtLeastNNotCoveredVertices here at all? we already have the shortest path, these computations seem redundant
                    newEdge = findShortestPath_AtLeastNNotCoveredVerticesMNewPixels(reebGraph, componentToVertices[selectedComponentId], terminals,
                                                                          vertexIsCoveredByRadius, N,
                                                                          fff,
                                                                          TotalPixelCoverage,
                                                                          MINIMUM_PIXEL_COVERAGE,
                                                                          origMask);
                    markVerticesCoveredByRadius(reebGraph, newEdge, componentToVertices, selectedComponentId, vertexIsCoveredByRadius);


                    if (!newEdge.empty()) {
                        tryToAddNewChain(reebGraph, origMask, newEdge, chains, TotalPixelCoverage);
                    }
                    else {
                        std::cout << "new edge is empty!\n";
                    }
                } while (!newEdge.empty());

                std::cout << "Chains: \n";
                for (int i = 0; i<chains.size(); i++)
                {
                    std::cout << "chain : " << i ;
                    printChain(chains[i], reebGraph, true);
                }

                std::cout << "\nRemoving non-valence-2 vertices \n -{\n";
//                solveJunctions(reebGraph, pts, selectedComponentId, chains, activePtsPerComponent);
//                glueAllJunctions(reebGraph, pts, selectedComponentId, chains, activePtsPerComponent);
                advancedGlueAllJunction(reebGraph, pts, selectedComponentId, chains, activePtsPerComponent);
                for (const auto &c : chains) {
                    //std::cout << " chain " << std::endl;
                    markVerticesCoveredByRadius(reebGraph, c, componentToVertices, selectedComponentId, vertexIsCoveredByRadius);
                }

                std::cout << "\n}-\n Done removing non-valence-2 vertices \n\n";


                for (const auto &c : chains) {
                    allowedChains.push_back(chainToVerticesSeq(c));
                    allowedChainToLabeledPointPair.push_back(
                            std::make_pair(
                                    getTheLabeledPointConnectedToMe(activePtsPerComponent, c.front().m_source, selectedComponentId),
                                    getTheLabeledPointConnectedToMe(activePtsPerComponent, c.back().m_target, selectedComponentId)
                            )
                    );
                }

                std::cout << "\n COMPONENT " << selectedComponentId << "\n End pts: ";
                for (auto v : activePtsPerComponent[selectedComponentId])
                    std::cout << "\n" << v.first << "\t (graph vertex " << v.second << "),\t coords: " << reebGraph[v.second].location.transpose();
                std::cout << std::endl;

                std::cout << "\n--- Selecting best curves! \n";
//            std::vector<std::vector<edge_descriptor>> newchains = selectBestCurves(reebGraph, chains, origMask, fff);
                std::vector<std::vector<edge_descriptor>> newchains = selectAllCurves(reebGraph, chains, origMask, fff);

                std::vector<MyPolyline> myPolys; // myPolys contains polylines in this component
                std::vector<MyPolyline> splitMyPolys;

                for (const auto& c : newchains)
                    myPolys.push_back(chainToPolyline(c, reebGraph));

                polylineresult.insert(polylineresult.end(), myPolys.begin(), myPolys.end());
                for (const auto& chain : newchains) {
                    finalchains.push_back(chainToVerticesSeq(chain));

                    for (vertex_descriptor v : {chain.front().m_source, chain.back().m_target}){
                        int labeledPointConnected = getTheLabeledPointConnectedToMe(activePtsPerComponent, v, selectedComponentId);
                        std::cout << "Vertex " << v << " is associated with labeled point " << labeledPointConnected << "; valence ++\n";
                        keypointToValenceMap[labeledPointConnected]++;
                    }
                }

            }
            else {
                std::cout << "The cluster turns out to be large. We should run a simpler procedure here!\n";
                // first we make a subgraph
                G subgraph = makeSubGraph(reebGraph, componentToVertices[selectedComponentId]);
                // second we break the subgraph to chains by high-valence vertices
                std::map<edge_descriptor, size_t> edge_to_simpleChain;

                std::map<edge_descriptor, size_t> ignore;
                removeBranchesFilter1(subgraph,false,ignore);

                std::vector<std::vector<typename G::edge_descriptor>> simpleChains = chainDecomposition(subgraph, edge_to_simpleChain);
                // third we add these chains and its ends to the output
                // upd allowedChainToLabeledPointPair
                // upd newlyAddedIntersections == newFarthestPointsToAdd_coords +
                // upd finalchains +
                // upd allowedChains +
                // upd keypointToValenceMap +
                // polylines
                std::vector<vertex_descriptor> simpleChainsEnds;
                for (const auto& chain : simpleChains) {
                    vertex_descriptor chain_start = chain.front().m_source;
                    vertex_descriptor chain_end = chain.back().m_target;
                    //check if current vertex is an active vertex for this component
                    auto chain_start_pointer = std::find_if(activePtsPerComponent[selectedComponentId].begin(),
                                                            activePtsPerComponent[selectedComponentId].end(),
                                                            [&chain_start](const std::pair<int, int>& el) {return el.second == chain_start; });
                    if (chain_start_pointer == activePtsPerComponent[selectedComponentId].end())
                    {
                        simpleChainsEnds.push_back(chain_start);
                        newFarthestPointsToAdd.push_back(chain_start);
                        newFarthestPointsToAdd_coords.emplace_back( reebGraph[chain_start].location.y(), reebGraph[chain_start].location.x() );
                        activePtsPerComponent[selectedComponentId].insert(std::make_pair(-newFarthestPointsToAdd.size(), chain_start));
                    }
                    auto chain_end_pointer = std::find_if(activePtsPerComponent[selectedComponentId].begin(),
                                                            activePtsPerComponent[selectedComponentId].end(),
                                                            [&chain_end](const std::pair<int, int>& el) {return el.second == chain_end; });
                    if (chain_start_pointer == activePtsPerComponent[selectedComponentId].end())
                    {
                        simpleChainsEnds.push_back(chain_end);
                        newFarthestPointsToAdd.push_back(chain_end);
                        newFarthestPointsToAdd_coords.emplace_back( reebGraph[chain_end].location.y(), reebGraph[chain_end].location.x() );
                        activePtsPerComponent[selectedComponentId].insert(std::make_pair(-newFarthestPointsToAdd.size(), chain_end));
                    }
                    finalchains.push_back(chainToVerticesSeq(chain));
                    allowedChains.push_back(chainToVerticesSeq(chain));
                    allowedChainToLabeledPointPair.push_back(
                            std::make_pair(
                                    getTheLabeledPointConnectedToMe(activePtsPerComponent, chain_start, selectedComponentId),
                                    getTheLabeledPointConnectedToMe(activePtsPerComponent, chain_end, selectedComponentId)
                            )
                    );
                    for (vertex_descriptor v : {chain.front().m_source, chain.back().m_target}){
                        int labeledPointConnected = getTheLabeledPointConnectedToMe(activePtsPerComponent, v, selectedComponentId);
                        std::cout << "Vertex " << v << " is associated with labeled point " << labeledPointConnected << "; valence ++\n";
                        keypointToValenceMap[labeledPointConnected]++;
                    }
                }

                std::vector<MyPolyline> myPolys;
                for (const auto& c : simpleChains)
                    myPolys.push_back(chainToPolyline(c, reebGraph));

                polylineresult.insert(polylineresult.end(), myPolys.begin(), myPolys.end());
                componentsToProcess.pop();
                continue;
            }
        }
        if (activePtsPerComponent[selectedComponentId].size() < 2){
//        if (activePtsPerComponent[selectedComponentId].size() == 1){
            assert(componentToVertices[selectedComponentId].size()>0);
            std::cout << "Component " << selectedComponentId << " has " << activePtsPerComponent[selectedComponentId].size() << " points. Trying to add some more\n";
            std::set<size_t> componentVertex;
            if (activePtsPerComponent[selectedComponentId].size() == 1) {
                // if we have 1 labeled point fill componentVertex with it
                for (const auto& a : activePtsPerComponent[selectedComponentId])
                    componentVertex.insert(a.second);
            }
            if (activePtsPerComponent[selectedComponentId].empty()) {
                // if no points labeled, fill it with first point we find
                componentVertex.insert(componentToVertices[selectedComponentId][0]);
                std::cout<<"0 actives, pushed a random "
                << componentToVertices[selectedComponentId][0] << " to find the farthest\n";
            }
            std::pair<bool, size_t> shouldIAddNewPointForThisComponent = findFarthestPixelUncoveredVertexFromTheSameComponent(reebGraph,
                                                                                                                              componentVertex,
                                                                                                                              TotalPixelCoverage,
                                                                                                                              MINIMUM_PIXEL_COVERAGE,
                                                                                                                              origMask);
            if (shouldIAddNewPointForThisComponent.first) {
                size_t newpoint_index = shouldIAddNewPointForThisComponent.second;
                // add this point to labeled points
                pts.emplace_back( reebGraph[newpoint_index].location.y(), reebGraph[newpoint_index].location.x() );
                newFarthestPointsToAdd_coords.emplace_back(reebGraph[newpoint_index].location.y(), reebGraph[newpoint_index].location.x());
                activePtsPerComponent[selectedComponentId].insert(std::make_pair(pts.size()-1, newpoint_index));
                std::cout << "\n ! Found new endpoint for component " << selectedComponentId << "\n";
                componentsToProcess.push(selectedComponentId);
            } else {
                std::cout << "\n ? Sorry, component " << selectedComponentId << " has " << activePtsPerComponent[selectedComponentId].size() << " active vertex but I can't find another one far enough \n";
            }
        }
        componentsToProcess.pop();
    } while (!componentsToProcess.empty());
    for (unsigned i = 0; i < polylineresult.size(); i++) {
        std::cout << i << " polyline: " << polylineresult[i][0].transpose() << " -> " << polylineresult[i].back().transpose() << std::endl;
    }
    newlyAddedIntersections = newFarthestPointsToAdd_coords;

    std::map<vertex_descriptor , int> endpoints_occurence;
    std::map<int, int> labelpoints_occurence;
    for (const auto& chain : finalchains) {
        endpoints_occurence[chain.front()]++;
        endpoints_occurence[chain.back()]++;
        labelpoints_occurence[getLabeledPointConnected(activePtsPerComponent, chain.front())]++;
        labelpoints_occurence[getLabeledPointConnected(activePtsPerComponent, chain.back())]++;
    }

    std::cout << "Labeled points:\n";
    for ( const auto &occ : labelpoints_occurence ) {
        std::cout << "(" << occ.first << ") : " << occ.second << "\n";
    }

    std::cout << "Active points:\n";
    for ( const auto &occ : endpoints_occurence ) {
        std::cout << "(" << getLabeledPointConnected(activePtsPerComponent, occ.first) << ") " << occ.first << " : " << occ.second << "\n";
    }

    std::cout<<" we didnt segfault YET "<<std::endl;
    return polylineresult;
}