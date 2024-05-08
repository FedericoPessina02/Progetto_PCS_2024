#pragma once
#include "Eigen/Eigen"
#include "vector"

using namespace std;

struct TracesMesh {
    vector<unsigned int> traces_id;
    vector<array<double,6>> traces_vertices;
    vector<array<unsigned int,2>> traces_fracture;
    vector<double> traces_length;
    TracesMesh();
    void addTrace();
};
