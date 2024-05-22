#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "classes/Fracture.hpp"
#include "classes/TracesMesh.hpp"

using namespace std;

namespace Algorithms {

map<int, vector<Fracture>> assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number);

void cutTracesInsidePartition(vector<Fracture>& fractures, TracesMesh& mesh);

void cutTracesOverlapping(vector<Fracture>& overlapping, vector<Fracture>& other_fractures, TracesMesh& mesh);

void cutTraces(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh, const int dimension);

void ordinaFract(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh, string nome_file);

}
