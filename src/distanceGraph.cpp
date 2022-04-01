#include "distanceGraph.h"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

std::vector<std::vector<double>> computeGraphDistancesPerPoly(const G& reebGraph, std::vector<MyPolyline> polys, const std::vector<std::vector<vertex_descriptor>>& allowedChains) {
    std::cout << "FIND Y JUNCTIONS" << std::endl;
    std::vector<size_t> interestingPoints; // points with high valence OR endpoints

    const int nbVert = boost::num_vertices(reebGraph);
    std::vector<size_t> valence(nbVert, 0);

    auto es = boost::edges(reebGraph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        //std::cout << "AN EDGE : "<< boost::source(*eit, reebGraph) << ' ' << boost::target(*eit, reebGraph) << std::endl;
        int source = boost::source(*eit, reebGraph);
        int target = boost::target(*eit, reebGraph);
        valence[source] += 1;
        valence[target] += 1;
    }

    for(int i=0; i<valence.size();i++)
    {
        if (valence[i] >= 3)
        {
            int countCurvesPassingBy = 0;
            for(int j=0; j<polys.size();j++)
            {
                for(int k=0; k<polys[j].size();k++)
                {
                    if (polys[j][k] == reebGraph[i].location)
                    {
                        //std::cout << "valence : " << valence[i] << " at pt : " << polys[j][k].x() << ", "
                        //          << polys[j][k].y() << ". Btw vertex index is : " << i << ", poly is " << j
                        //          << ", and poly vert index is " << k << std::endl;
                        countCurvesPassingBy++;
                        break;
                    }
                }
            }
            //std::cout<<"nb of curves passing by vertex "<<i<<" : "<<countCurvesPassingBy<<std::endl;
            if (countCurvesPassingBy >=2)
            {
                interestingPoints.push_back(i);
            }
        }
    }

    for (const auto& c : allowedChains) {
        interestingPoints.push_back(c.front());
        interestingPoints.push_back(c.back());
    }

//    for(int i=0; i<valence.size();i++)
//    {
//        if (valence[i] >= 3)
//        {
//            interestingPoints.push_back(i);
//        }
//    }

    //std::cout<<" Final list of y junctions."<<std::endl;
    for (auto e : interestingPoints)
    {
        //std::cout<< "vertex "<<e<<" at position "<<reebGraph[e].location.x()<< ", "<<reebGraph[e].location.y()<<std::endl;
    }

    if (interestingPoints.size() == 0)
    {
        std::vector<std::vector<double>> minDistPerPolyPerVertex;
        for (auto poly : polys)
        {
            std::vector<double> tmp_distances;
            for (auto v : poly)
            {
                tmp_distances.push_back(0.1);
            }
            minDistPerPolyPerVertex.push_back(tmp_distances);
        }
        return minDistPerPolyPerVertex;
    }

    std::vector<vertex_descriptor> pDij(num_vertices(reebGraph));
    std::vector<double> mindDij(num_vertices(reebGraph));
    auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, reebGraph));
    auto distMap = make_iterator_property_map(mindDij.begin(), get(&Cluster::clusterIdx, reebGraph));
    boost::dijkstra_shortest_paths(reebGraph, interestingPoints[0],predecessor_map(predMap).distance_map(distMap).weight_map(get(&Edge::weight, reebGraph)));
    //std::cout<<"y junctions size"<<interestingPoints.size()<<std::endl;

    //for(int y=0; y<num_vertices(reebGraph); y++)
    for (int yIdx = 1; yIdx < interestingPoints.size(); yIdx ++)
    {
        //auto y = interestingPoints[yIdx];
        //std::cout<<" y for dist"<<interestingPoints[yIdx]<<std::endl;
        // Do shortest path to that y junction
        std::vector<vertex_descriptor> pDij(num_vertices(reebGraph));
        std::vector<double> dDij(num_vertices(reebGraph));
        auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, reebGraph));
        auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx, reebGraph));
        boost::dijkstra_shortest_paths(reebGraph, interestingPoints[yIdx],predecessor_map(predMap).distance_map(distMap).weight_map(get(&Edge::weight, reebGraph)));
        //std::cout<<dDij[interestingPoints[yIdx]]<<std::endl;
//        if (y == 1218)
//        {
//            for (int i = 0; i<dDij.size(); i++) {
//                //std::cout<<"vert "<<i<< " coords : "<<reebGraph[i].location.x()<< ", "<<reebGraph[i].location.y()<<" dist : "<<dDij[i]<<std::endl;
//            }
//        }


        for (int i =0; i<dDij.size(); i++)
        {
            dDij[i] < mindDij[i]? mindDij[i] = dDij[i]: mindDij[i] = mindDij[i];
        }
    }
    for (int i : {18,19,20,21})
    {
        //std::cout<< "poly nb "<<i<<std::endl;
        for (auto e : polys[i])
        {
            //std::cout<<e.x()<<" "<<e.y()<<std::endl;
        }
    }
    std::vector<std::vector<double>> minDistPerPolyPerVertex;
    for (auto poly : polys)
    {
        std::vector<double> tmp_distances;
        for (auto v : poly)
        {
            bool foundAMatch = false;
            double minDistVert;
            for(int i = 0 ; i<mindDij.size(); i++)
            {
                if (!foundAMatch && reebGraph[i].location == v)
                {
                    minDistVert = mindDij[i];
                    foundAMatch = true;
                }
                else if (reebGraph[i].location == v && mindDij[i] < minDistVert) minDistVert = mindDij[i];
            }
            if (!foundAMatch) std::cout<< " HAVENT FOUND A MATCH FOR VERTEX "<< v <<std::endl;
            else
            {
                tmp_distances.push_back(minDistVert);
            }

            //std::cout<<"mindisties "<<minDistVert<<std::endl;
        }
        minDistPerPolyPerVertex.push_back(tmp_distances);
    }
    for (int i = 0 ;i<polys.size();i++)
    {
        //std::cout<<"POLY "<< i<<std::endl;
        for (int j =0; j<polys[i].size(); j++)
        {
            //std::cout<<polys[i][j].x()<<" "<<polys[i][j].y()<< " dist : "<<minDistPerPolyPerVertex[i][j]<<std::endl;
        }
    }



    return minDistPerPolyPerVertex;
}
