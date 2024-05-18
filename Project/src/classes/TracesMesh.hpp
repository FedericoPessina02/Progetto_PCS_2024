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
    void addTrace(unsigned int& id, array<Eigen::Vector3d,2>& vertices, array<unsigned int,2>& fractures) {
        traces_id.push_back(id);
        traces_vertices.push_back(vertices);
        traces_fracture.push_back(fractures);
        traces_length.push_back((vertices[0] - vertices[1]).norm());
    }
};
