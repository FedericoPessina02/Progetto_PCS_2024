#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "classes/Fracture.hpp"
#include "classes/TracesMesh.hpp"

using namespace std;

namespace Algorithms {

vector<Eigen::Vector3d> calculateIntersectionsPoints(Fracture& fracture, Eigen::Vector3d line, Eigen::Vector3d point);

map<int, vector<Fracture>> assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number);

void cutTracesInsidePartition(vector<Fracture>& fractures, TracesMesh& mesh);

void cutTracesOverlapping(vector<Fracture>& overlapping, vector<Fracture>& other_fractures, TracesMesh& mesh);

void cutTraces(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh);

vector<PolygonalMesh> cutPolygonalMesh(map<int, vector<Fracture>>& id_to_fractures);

void cutPolygonBySegment(PolygonalMesh& mesh, unsigned int polygonId, array<Eigen::Vector3d,2> segment);

}
