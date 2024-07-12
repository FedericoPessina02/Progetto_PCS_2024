#pragma once
#include <iostream>
#include <vector>
#include "classes/Fracture.hpp"
#include "Eigen/Eigen"
#include "PolygonalMesh.hpp"

using namespace std;

namespace Utils{

inline int tol_coeff = 500;

vector<Fracture> fractureInput(const string& filename, array<double, 6>& domain_borders);

vector<Eigen::Vector3d> calculateDistinctPoints(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

bool compareSegments(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b);

// le performance non sembrano migliorare quindi non l'abbiamo inserito

// inline vector<Eigen::Vector3d> calculateDistinctPoints(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b) {
//     // la funzione verifica che a e b siano due coppie di punti distinti
//     vector<Eigen::Vector3d> result;
//     result.push_back(a[0]);
//     result.push_back(a[1]);
//     for (const Eigen::Vector3d& b_el : b) {
//         if ((a[0]-b_el).norm() >= tol_coeff*numeric_limits<double>::epsilon() && (a[1]-b_el).norm() >= tol_coeff*numeric_limits<double>::epsilon()) {
//             result.push_back(b_el);
//         }
//     }
//     return result;
// }

// inline bool compareSegments(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b) {
//     // verifico se i due segmenti sono in realtà lo stesso a meno di permutazione degli elementi
//     if ((a[0]-b[0]).squaredNorm() < tol_coeff*numeric_limits<double>::epsilon() && (a[1]-b[1]).squaredNorm() < tol_coeff*numeric_limits<double>::epsilon()) {
//         return true;
//     }
//     if ((a[0]-b[1]).squaredNorm() < tol_coeff*numeric_limits<double>::epsilon() && (a[1]-b[0]).squaredNorm() < tol_coeff*numeric_limits<double>::epsilon()) {
//         return true;
//     }
//     return false;
// }


void Stampa1(string nome_file,TracesMesh& mesh);

void ExportSTL(string nome_file, vector<PolygonalMesh>& mesh_collection);

void Merge(vector<pair<int, double>>& v, const unsigned int& sx, const unsigned int& cx, const unsigned int& dx);

void MergeSort(vector<pair<int, double>>& v);

void MergeSort(vector<pair<int, double>>& v, const unsigned int& sx, const unsigned int& dx);

void BubbleSort(std::vector<pair<int, double>>& data);

}
