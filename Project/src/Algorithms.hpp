#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "classes/Fracture.hpp"
#include "classes/TracesMesh.hpp"

using namespace std;

namespace Algorithms {

void assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number);

void cutTraces(vector<Fracture>& fractures, TracesMesh& mesh);

}
