#pragma once
#include "Eigen/Eigen"
#include "vector"

using namespace std;

struct TracesMesh {
    vector<unsigned int> traces_id;
    vector<array<Eigen::Vector3d,2>> traces_vertices;
    vector<array<unsigned int,2>> traces_fracture;
    vector<double> traces_length;
    TracesMesh() = default;
    void addTrace(unsigned int& id, array<double,6>& vertices, vector<array<unsigned int,2>>& traces_fracture);
};
