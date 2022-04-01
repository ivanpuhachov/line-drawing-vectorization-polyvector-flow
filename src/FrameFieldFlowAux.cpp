//
// Created by william on 2020-12-11.
//

#include "stdafx.h"
#include "FrameFieldFlow.h"
#include <complex>
#include <Eigen/Sparse>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323
#endif 

void analysis(const MyPolyline& poly)
{
    std::cout << "Starting polyline analysis" << std::endl;
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    auto max_ele = std::max_element(ll.begin(), ll.end());
    std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    std::sort(ll.begin(), ll.end());

    for (auto i : poly)
    {
        //    std::cout << i[0] << " ; " << i[1] << std::endl;
    }

    for (auto i : ll)
    {
        //    std::cout << i << std::endl;
    }

    int lN = ll.size();
    double total_length = accumulate(ll.begin(), ll.end(), 0.0);
    double average = total_length/lN;

    std::cout << "length of ll (list of edge length) : " << lN << std::endl;
    if (true)
    {
        std::cout << "min length : " << ll[0] << std::endl;
        std::cout << "10th percentile : " << ll[int(lN / 10)] << std::endl;
        std::cout << "20th percentile : " << ll[int(2 * lN / 10)] << std::endl;
        std::cout << "30th percentile : " << ll[int(3 * lN / 10)] << std::endl;
        std::cout << "40th percentile : " << ll[int(4 * lN / 10)] << std::endl;
        std::cout << "50th percentile : " << ll[int(5 * lN / 10)] << std::endl;
        std::cout << "60th percentile : " << ll[int(6 * lN / 10)] << std::endl;
        std::cout << "70th percentile : " << ll[int(7 * lN / 10)] << std::endl;
        std::cout << "80th percentile : " << ll[int(8 * lN / 10)] << std::endl;
        std::cout << "90th percentile : " << ll[int(9 * lN / 10)] << std::endl;
    }
    std::cout << "max length : " << ll[int(lN - 1)] << std::endl;
    std::cout << "Average edge length : " << average << std::endl;
    std::cout << "Total curve length : " << total_length << std::endl;

    std::cout << "Ending polyline analysis\n" << std::endl;

}

double total_edge_length(const MyPolyline& poly)
{
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    double total_length = accumulate(ll.begin(), ll.end(), 0.0);
    return total_length;
}

int min_edge_index(const MyPolyline& poly)
{
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    auto min_ele = std::min_element(ll.begin(), ll.end());
    //std::cout<<" MIN ELEMENT INDEX : " << min_ele - ll.begin()<<std::endl;
    return min_ele - ll.begin();
}

double min_edge_length(const MyPolyline& poly)
{
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    auto min_ele = std::min_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT : " << ll[min_ele - ll.begin()]<<std::endl;
    return ll[min_ele - ll.begin()];
}

int max_edge_index(const MyPolyline& poly)
{
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    auto max_ele = std::max_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    return max_ele - ll.begin();
}

double max_edge_length(const MyPolyline& poly)
{
    int N = poly.size();

    // edges length
    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }
    auto max_ele = std::max_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    return ll[max_ele - ll.begin()];
}

MyPolyline reverse(const MyPolyline& poly)
{
    std::vector<Eigen::Vector2d> newPoly;
    for (int i = poly.size()-1; i >= 0; i--)
    {
        newPoly.push_back(poly[i]);
    }
    return newPoly;
}

std::vector<double> angles(const MyPolyline& poly)
{
    int N = poly.size();

    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }

    std::vector<double> a;
    for (int i = 1; i < N - 1; i++)
    {
        double ok = ((poly[i][0] - poly[i-1][0]) * (poly[i+1][0] - poly[i][0]) + (poly[i][1] - poly[i-1][1]) * (poly[i+1][1] - poly[i][1])) / (ll[i-1] * ll[i]);
        if (ok >= 1)
        {
            //std::cout<<i<<std::endl;
            a.push_back(0.0);
        }
        else {
            double ang = acos(
                    ((poly[i][0] - poly[i - 1][0]) * (poly[i + 1][0] - poly[i][0]) +
                     (poly[i][1] - poly[i - 1][1]) * (poly[i + 1][1] - poly[i][1])) / (ll[i - 1] * ll[i])
            );
            //std::cout << "length i-1 : " << ll[i - 1] << " , length i : " << ll[i] << std::endl;
            //std::cout << "a x : " << poly[i][0] - poly[i - 1][0] << " , a y : " << poly[i][1] - poly[i - 1][1]<< std::endl;
            //std::cout << "b x : " << poly[i + 1][0] - poly[i][0] << " , b y : " << poly[i + 1][1] - poly[i][1]<< std::endl;
            //std::cout << "thing to be acosed: " << ((poly[i][0] - poly[i - 1][0]) * (poly[i + 1][0] - poly[i][0]) +
              //                                      (poly[i][1] - poly[i - 1][1]) * (poly[i + 1][1] - poly[i][1])) /
                //                                   ll[i - 1] / ll[i] << std::endl;
            //std::cout << "angle ? : " << ang * 360 / 2 / M_PI << std::endl;
            a.push_back(ang);
        }
    }
    //for (int i =0; i<a.size(); i++){
     //   std::cout<<"ANGLE "<<i<<" : "<<a[i]<<std::endl;
    //}
    int lN = a.size();
    double total_angle = accumulate(a.begin(), a.end(), 0.0);
    double average = total_angle/lN;

    std::sort(a.begin(), a.end());

    //std::cout << "length of a (list of edge length) : " << lN << std::endl;
    if (false)
    {
        std::cout << "min angle : " << a[0] << std::endl;
        std::cout << "10th percentile : " << a[int(lN / 10)] << std::endl;
        std::cout << "20th percentile : " << a[int(2 * lN / 10)] << std::endl;
        std::cout << "30th percentile : " << a[int(3 * lN / 10)] << std::endl;
        std::cout << "40th percentile : " << a[int(4 * lN / 10)] << std::endl;
        std::cout << "50th percentile : " << a[int(5 * lN / 10)] << std::endl;
        std::cout << "60th percentile : " << a[int(6 * lN / 10)] << std::endl;
        std::cout << "70th percentile : " << a[int(7 * lN / 10)] << std::endl;
        std::cout << "80th percentile : " << a[int(8 * lN / 10)] << std::endl;
        std::cout << "90th percentile : " << a[int(9 * lN / 10)] << std::endl;
    }
    std::cout << "max angle : " << a[int(lN - 1)] * 360 / 2/M_PI<< std::endl;
    std::cout << "Average edge angle : " << average * 360 / 2/M_PI<< std::endl;
    std::cout << "Total curve angle : " << total_angle * 360 / 2/M_PI<< std::endl;

    //auto max_ele = std::max_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    //return ll[max_ele - ll.begin()];


    //std::cout << "Ending polyline analysis\n" << std::endl;
    return a;
}

double max_angle(const MyPolyline& poly)
{
    int N = poly.size();

    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }

    std::vector<double> a;
    for (int i = 1; i < N - 1; i++)
    {
        double ok = ((poly[i][0] - poly[i-1][0]) * (poly[i+1][0] - poly[i][0]) + (poly[i][1] - poly[i-1][1]) * (poly[i+1][1] - poly[i][1])) / (ll[i-1] * ll[i]);
        if (ok >= 1)
        {
            //std::cout<<i<<std::endl;
            a.push_back(0.0);
        }
        else {
            double ang = acos(
                    ((poly[i][0] - poly[i - 1][0]) * (poly[i + 1][0] - poly[i][0]) +
                     (poly[i][1] - poly[i - 1][1]) * (poly[i + 1][1] - poly[i][1])) / (ll[i - 1] * ll[i])
            );
            a.push_back(ang);
        }
    }

    int lN = a.size();

    std::sort(a.begin(), a.end());
    double max_angle =  a[int(lN - 1)] * 360 / 2/M_PI;
    //std::cout << "max angle : " <<max_angle<< std::endl;

    //auto max_ele = std::max_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    //return ll[max_ele - ll.begin()];


    //std::cout << "Ending polyline analysis\n" << std::endl;
    return max_angle;
}

double neighbour_avg_angle(const MyPolyline& poly, int k)
{
    int N = poly.size();


    std::vector<double> ll = {};
    for (int i = 1; i < N; i++)
    {
        ll.push_back(sqrt(pow(poly[i-1][0] - poly[i][0], 2) + pow(poly[i-1][1] - poly[i][1], 2)));
    }

    std::vector<double> a;
    for (int i = 1; i < N - 1; i++)
    {
        double ok = ((poly[i][0] - poly[i-1][0]) * (poly[i+1][0] - poly[i][0]) + (poly[i][1] - poly[i-1][1]) * (poly[i+1][1] - poly[i][1])) / (ll[i-1] * ll[i]);
        if (ok >= 1)
        {
            a.push_back(0.0);
        }
        else {
            double ang = acos(
                    ((poly[i][0] - poly[i - 1][0]) * (poly[i + 1][0] - poly[i][0]) +
                     (poly[i][1] - poly[i - 1][1]) * (poly[i + 1][1] - poly[i][1])) / (ll[i - 1] * ll[i])
            );
            a.push_back(ang);
        }
    }
    int lN = a.size();
    double total_angle = accumulate(a.begin(), a.end(), 0.0);
    double average = total_angle/lN;

    //std::sort(a.begin(), a.end());

    //std::cout << "max angle : " << a[int(lN - 1)] * 360 / 2/M_PI<< std::endl;
    //std::cout << "Average edge angle : " << average * 360 / 2/M_PI<< std::endl;
    //std::cout << "Total curve angle : " << total_angle * 360 / 2/M_PI<< std::endl;

    //auto max_ele = std::max_element(ll.begin(), ll.end());
    //std::cout<<" MAX ELEMENT INDEX : " << max_ele - ll.begin()<<std::endl;
    //return ll[max_ele - ll.begin()];


    //std::cout << "Ending polyline analysis\n" << std::endl;

    auto max_ele_pointer = std::max_element(a.begin(), a.end());
    auto max_ele_idx = max_ele_pointer - a.begin();
    auto max_ele = a[max_ele_idx];

    int min_range_idx = max_ele_idx - k >=0 ? max_ele_idx - k : 0;
    int max_range_idx = max_ele_idx + k <= a.size() - 1  ? max_ele_idx + k : a.size() - 1;

    double angle_acc = 0.0;
    for(int i = min_range_idx; i<= max_range_idx; ++i)
    {
        angle_acc += a[i];
    }
    angle_acc = angle_acc / (max_range_idx - min_range_idx + 1);

    // Maybe i should return max_range_idx -  min_range_idx + 1, so that i know how many angles went into the avg.

    return angle_acc;

    return 0.0;
}