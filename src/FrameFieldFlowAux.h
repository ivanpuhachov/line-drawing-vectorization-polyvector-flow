//
// Created by william on 2020-12-11.
//

#ifndef POLYVECTOR_THING_FRAMEFIELDFLOWAUX_H
#define POLYVECTOR_THING_FRAMEFIELDFLOWAUX_H

#pragma once
#include "typedefs.h"


double total_edge_length(const MyPolyline& poly);
int min_edge_index(const MyPolyline& poly);
double min_edge_length(const MyPolyline& poly);
int max_edge_index(const MyPolyline& poly);
double max_edge_length(const MyPolyline& poly);
void analysis(const MyPolyline& poly);
MyPolyline reverse(const MyPolyline& poly);
std::vector<double> angles(const MyPolyline& poly);
double max_angle(const MyPolyline& poly);
double neighbour_avg_angle(const MyPolyline& poly, int k);


#endif //POLYVECTOR_THING_FRAMEFIELDFLOWAUX_H
