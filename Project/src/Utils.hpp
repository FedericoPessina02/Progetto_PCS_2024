#pragma once
#include <iostream>
#include <vector>
#include "classes/Fracture.hpp"
#include "Eigen/Eigen"

using namespace std;

namespace Utils{

inline int tol_coeff = 500;

vector<Fracture> fractureInput(const string& filename, array<double, 6>& domain_borders);

vector<Eigen::Vector3d> calculateDistinctPoints(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

bool compareSegments(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

void Stampa1(string nome_file,TracesMesh& mesh);

void ExportSTL(string nome_file, vector<PolygonalMesh>& mesh_collection);

void Merge(vector<pair<int, double>>& v, const unsigned int& sx, const unsigned int& cx, const unsigned int& dx);

void MergeSort(vector<pair<int, double>>& v);

void MergeSort(vector<pair<int, double>>& v, const unsigned int& sx, const unsigned int& dx);

void BubbleSort(std::vector<pair<int, double>>& data);

}
