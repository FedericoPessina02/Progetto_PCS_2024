#pragma once
#include <iostream>
#include <vector>
#include "classes/Fracture.hpp"
#include "Eigen/Eigen"

using namespace std;

namespace Utils{

vector<Fracture> fractureInput(const string& filename, array<double, 6>& domain_borders);

vector<Eigen::Vector3d> calculateDistinctPoints(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

bool compareSegments(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

void Stampa1(string nome_file,TracesMesh& mesh);

void ExportSTL(string nome_file, vector<PolygonalMesh>& mesh_collection);

}
