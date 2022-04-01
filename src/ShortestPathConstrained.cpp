#include "ShortestPathConstrained.h"

bool yJunctionTest(const std::vector<vertex_descriptor>& curPath, vertex_descriptor newV, const G& g)
{
    // curPath.back() is the intersection vertex.
    // curPath.front() is a neighbor from the 'principal' branch
    // newV is a neighbor from 'another' branch.
    if ((boost::degree(curPath.back(), g) < 3) || curPath.size()<2 /*|| g[curPath.back()].clusterCurveHitSingularity || g[curPath.back()].nextToSingularity*/ || (g[curPath.back()].seedCurve==-1)) //seedCurve = -1 for the added junctions
        return true;
    //find all shared traced curves
    //std::set<size_t> sharedTracedCurves;
    std::map<size_t, std::set<size_t>> myCurves;
    std::vector<size_t> testTriplet = { curPath[curPath.size() - 2], curPath.back(), newV }; // WN : is curPath ever bigger than 2 vertices? Putting all 3 vertices in the same variable, intersection vertex in the middle

    //std::cout<<"hello from y junction test ! "<<testTriplet[0]<<" "<<testTriplet[1]<<" "<<testTriplet[2]<<std::endl;

    for (size_t v : testTriplet) // WN : iterate over the 3 vertices.
    {
        for (const auto& p : g[v].clusterPoints) // WN : p is a const PointOnCurve ? ***
            myCurves[v].insert(p.curve); // WN: What is p and p.curve ? ***

    }

    std::set<size_t> sharedTracedCurves01, sharedTracedCurves12; // WN : Intersection between what is in myCurves[v : testTriplet]
    std::set_intersection(myCurves[testTriplet[0]].begin(),
                          myCurves[testTriplet[0]].end(),
                          myCurves[testTriplet[1]].begin(),
                          myCurves[testTriplet[1]].end(),
                          std::inserter(sharedTracedCurves01,
                                        sharedTracedCurves01.end()));
    std::set_intersection(myCurves[testTriplet[1]].begin(),
                          myCurves[testTriplet[1]].end(),
                          myCurves[testTriplet[2]].begin(),
                          myCurves[testTriplet[2]].end(),
                          std::inserter(sharedTracedCurves12,
                                        sharedTracedCurves12.end()));

    //now "orient" the valence-3 cluster curve
    std::map<size_t, std::array<PointOnCurve, 2>> ptsOnSharedCurves01, ptsOnSharedCurves12;
    for (int idx=0; idx<2; ++idx)
    {
        size_t v = testTriplet[idx];
        for (const auto& p : g[v].clusterPoints)
        {
            if (sharedTracedCurves01.find(p.curve) != sharedTracedCurves01.end()) // WN : If p.curve (wtv it is..) is part of sharedTracedCurves01 (which is a bunch on indices, of what?), add p to
                ptsOnSharedCurves01[p.curve][idx] = p;
        }
    }

    Eigen::Vector2d root = g[curPath.back()].root;
    int shouldFlipRoot = 0;
    for (const auto& vv : ptsOnSharedCurves01)
    {
        Eigen::Vector2d vec = vv.second[1].p - vv.second[0].p;
        if (vec.dot(root) < 0)
            shouldFlipRoot++;
        else
            shouldFlipRoot--;
    }
    if (shouldFlipRoot > 0)
        root = -root;

    for (int idx = 1; idx < 3; ++idx)
    {
        size_t v = testTriplet[idx];
        for (const auto& p : g[v].clusterPoints)
        {
            if (sharedTracedCurves12.find(p.curve) != sharedTracedCurves12.end())
                ptsOnSharedCurves12[p.curve][idx-1] = p;
        }
    }

    int isOrientationRight = 0;
    for (const auto& vv : ptsOnSharedCurves12)
    {
        Eigen::Vector2d vec = vv.second[1].p - vv.second[0].p;
        if (vec.dot(root) < 0)
            isOrientationRight--;
        else
            isOrientationRight++;
    }

    //std::cout << "yay";
    return isOrientationRight >= 0;
}

size_t countSharedCurves(vertex_descriptor vertex1, vertex_descriptor vertex2, const G& g){
    std::map<size_t, std::set<size_t>> myCurves;
    for (const auto& p : g[vertex1].clusterPoints)
        myCurves[vertex1].insert(p.curve);
    for (const auto& p : g[vertex2].clusterPoints)
        myCurves[vertex2].insert(p.curve);
    std::set<size_t> sharedTracedCurves; // WN : Intersection between what is in myCurves[v : testTriplet]
    std::set_intersection(myCurves[vertex1].begin(),
                          myCurves[vertex1].end(),
                          myCurves[vertex2].begin(),
                          myCurves[vertex2].end(),
                          std::inserter(sharedTracedCurves,
                                        sharedTracedCurves.end()));
    return sharedTracedCurves.size();
}

std::vector<edge_descriptor> findShortestPath_AtLeastNNotCoveredVerticesMNewPixels(
        const G& g,
        const std::vector<vertex_descriptor>& myComponent,
        const std::set<vertex_descriptor>& terminals,
        const std::map<vertex_descriptor,bool>& covered,
        const int N,
        const FrameFieldFlow& fff,
        const std::set<std::pair<int, int>>& existingCoverage,
        const int M,
        const cv::Mat& mask)
{
    //algorithm from https://link.springer.com/article/10.1007/s10707-019-00385-8
    struct Label
    {
        int numberOfUncoveredVerts;
        double cost;
        std::vector<vertex_descriptor> path;
    };

    typedef std::pair<vertex_descriptor,Label> VertexWithLabel;
    auto cmp = [](const VertexWithLabel& left, const VertexWithLabel& right) {return left.second.cost < right.second.cost; };

    auto isUncovered = [&covered](vertex_descriptor v)
    {
        auto it = covered.find(v);
        return (it == covered.end() || !it->second);
    };

    //std::ofstream f("output1.txt");
    std::map<std::pair<vertex_descriptor,vertex_descriptor>, std::vector<Label>> pathsBetweenTerminals;
    if (terminals.size() > 200) // TODO: optimize the code so it works with big clusters
    {
        std::cout << "WARNING: This code has optimization issues when cluster has" << terminals.size() << " (too many) terminals. Exiting!\n";
        return {};
    }
    for (vertex_descriptor t : terminals)
    {
        std::vector<VertexWithLabel> Q;
        //std::vector< VertexWithLabel>& qVec = Container<VertexWithLabel, std::vector<VertexWithLabel>, decltype(cmp)>(Q);

        std::map<vertex_descriptor, Label> bestLabels;
        for (auto v : myComponent)
        {
            bestLabels[v].cost = std::numeric_limits<double>::max();
            bestLabels[v].numberOfUncoveredVerts = (isUncovered(v) ? 1 : 0);
            bestLabels[v].path = {};
        }

        bestLabels[t].cost = 0;
        bestLabels[t].numberOfUncoveredVerts = 0; //we know the terminals are covered
        bestLabels[t].path = { t };


        std::map<vertex_descriptor,std::vector<Label>> Out_Lbls; //for each vertex and number of uncovered verts, we know the min dist
        //now find shortest paths from t
        Q = { { t,bestLabels[t]} };
        while (!Q.empty())
        {
            const auto& vIt = std::min_element(Q.begin(), Q.end(),cmp);

            vertex_descriptor u = vIt->first;
            auto curPath = vIt->second.path;
            double curPathLength = vIt->second.cost;
            int curNumOfUncoveredVerts = vIt->second.numberOfUncoveredVerts;
            
            Out_Lbls[u].push_back(vIt->second);
            
            Q.erase(vIt);

            for (auto [eit, eend] = boost::out_edges(u, g); eit != eend; ++eit)
            {
                vertex_descriptor v = eit->m_target;
                if (std::find(curPath.begin(), curPath.end(), v) != curPath.end()) //we only want simple paths
                {
                    continue;
                }
                if (g[v].nextToSingularity &&!yJunctionTest(curPath, v, g)){ // when vertex is near singularity yJunctionTest is unstable
                    continue;
                }

                Eigen::Vector2d yxSource(g[eit->m_source].location.y(), g[eit->m_source].location.x()), yxTarget(g[eit->m_target].location.y(), g[eit->m_target].location.x());
                double newLength = curPathLength + g[*eit].weight;//fff.energyPerEdge(yxSource, yxTarget); 
                int newNumOfUncoveredVerts = curNumOfUncoveredVerts;
                auto it = covered.find(v);
                if (it == covered.end() || !it->second)
                    newNumOfUncoveredVerts++;

                Label newLabel;
                newLabel.cost = newLength;
                newLabel.numberOfUncoveredVerts = newNumOfUncoveredVerts;
                newLabel.path = curPath;
                newLabel.path.push_back(v);

                if ((newLength < bestLabels[v].cost) && (newNumOfUncoveredVerts >= bestLabels[v].numberOfUncoveredVerts))
                {
                    //our new lengths/# are strictly better than what's stored, update the record
                    bestLabels[v] = newLabel;

                    bool alreadyInQueue = false;
                    for (auto& it: Q)
                    {
                        if (it.first == v)
                        {
                            it.second = newLabel; //update the existing label
                            alreadyInQueue = true;
                        }
                    }
                    if (!alreadyInQueue)
                        Q.push_back({ v, newLabel });
                }
                else if (((newLength+1e-6 < bestLabels[v].cost) && (newNumOfUncoveredVerts < bestLabels[v].numberOfUncoveredVerts)) || ((newLength > bestLabels[v].cost) && (newNumOfUncoveredVerts > bestLabels[v].numberOfUncoveredVerts)))
                {
                    //mixed case, have to duplicate the vertex
                    //look for a dominated vertex in the queue

                    //do I need that??
                    bool foundDominated = false;
                    for (const auto& vIt : Q)
                    {
                        if ((vIt.first == u) && ((vIt.second.cost < newLength + 1e-6) || (vIt.second.numberOfUncoveredVerts >= newNumOfUncoveredVerts)))
						{
							foundDominated = true;
							break;
						}
                    }

                    if (!foundDominated)
                    {
                        //push the duplicate
                        Q.push_back({ v,newLabel });
                    }
                }
            }
        }


        for (auto t1 : terminals)
        {
            //if (t1 < t)
            //{
            //TODO: somewhere here the is a bug: even though the graph is unoriented, somehow paths{t, t1} != paths{t1, t} 
                pathsBetweenTerminals[{t, t1}] = Out_Lbls[t1];
            //}
        }
        
        /*f    << "For the terminal " << t << " I found the following paths:" << std::endl;
        for (const auto& it : Out_Lbls)
        {
            f << "Vertex " << it.first << std::endl;
            for (const auto& it2 : it.second)
            {
                f << it2.cost << ", " << it2.numberOfUncoveredVerts << " ";
                for (auto v : it2.path)
                    f << v << "  ";
                f << std::endl;
            }
            f << std::endl;
        }*/
        //break;
    }

    //now choose the min length path connecting any of the terminals
    double bestPathCost = 1e10;
    int bestPixelCoverage = 0;
    std::vector<vertex_descriptor> bestPath;
    for (const auto& records : pathsBetweenTerminals)
    {
        for (const auto& label : records.second)
        {
            if ((label.cost < bestPathCost) && (label.numberOfUncoveredVerts > N))
            {
                int newCoverage = computeAddedCoverage(g, existingCoverage, label.path, mask);
                if ((newCoverage > bestPixelCoverage) && (newCoverage >= M)) {
                    bestPath = label.path;
                    bestPathCost = label.cost;
                    bestPixelCoverage = newCoverage;
                }
            }
        }
    }

    if (bestPath.size() > 0)
    {
        std::cout << "BEST PATH: " << bestPathCost << std::endl;
        std::cout << "BEST PATH PIXEL COVERAGE " << bestPixelCoverage << std::endl;
    }

    //FOR DEBUG
    auto isInside = [](const std::vector<size_t>& vec, size_t val)
    {
        return std::find(vec.begin(), vec.end(), val) != vec.end();
    };

    //std::vector<size_t> interesting = {1517, 492, 183};
    for (auto t: terminals)
    {
        for (auto t1 : terminals)
        {
            if (t1 < t)
            {
                /*if (isInside(interesting, t) && isInside(interesting, t1))
                {
                    std::cout << "Shortest path between " << t << " and " << t1 << ": ";
                    for (const auto& label : pathsBetweenTerminals[{t, t1}])
                    {
                        std::cout << "COST= " << label.cost << ", #uncovered vertices: " << label.numberOfUncoveredVerts << std::endl;
                        std::cout << "PATH= ";
                        for (int i = 0; i + 1 < label.path.size(); ++i)
                        {
                            Eigen::Vector2d yxSource(g[label.path[i]].location.y(), g[label.path[i]].location.x()), yxTarget(g[label.path[i + 1]].location.y(), g[label.path[i + 1]].location.x());
                            std::cout << label.path[i] << " (" << fff.energyPerEdge(yxSource, yxTarget) << ") - ";
                        }

                    }
                    std::cout << std::endl;
                }*/
            }
        }
    }

    std::vector<edge_descriptor> result;
    if (bestPath.size() > 1)
    {
        std::cout << std::endl;
        std::cout << "Best path cost: " << bestPathCost << std::endl;
        for (int i = 0; i+1 < bestPath.size(); ++i)
        {
            auto e = boost::edge(bestPath[i], bestPath[i + 1],g);
            result.push_back(e.first);
        }
    }
    return result;
}
