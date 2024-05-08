#pragma once
#include "Eigen/Eigen"
#include "vector"

using namespace std;

struct PolygonalMesh {
    unsigned int FractureId;
    vector<unsigned int> IdCell0Ds;
    vector<Eigen::Vector3d> CoordinateCell0Ds;
    vector<unsigned int> IdCell1Ds;
    vector<array<unsigned int,2>> VerticesCell1Ds;
    vector<unsigned int> IdCell2Ds;
    vector<vector<unsigned int>> VerticesCell2Ds;
    vector<vector<unsigned int>> EdgesCell2Ds;
    PolygonalMesh() = default;
};
