#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "classes/Fracture.hpp"
#include "classes/TracesMesh.hpp"

using namespace std;

namespace Algorithms {

void assignPartition(vector<Fracture>& fractures, map<int, vector<Fracture>>& id_to_fractures, array<double, 6>& domain_borders, const int partitions_number);

void cutTracesInsidePartition(vector<Fracture>& fractures, TracesMesh& mesh);

void cutTracesOverlapping(vector<Fracture>& overlapping, vector<Fracture>& other_fractures, TracesMesh& mesh);

void cutTraces(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh, const int dimension);

}
