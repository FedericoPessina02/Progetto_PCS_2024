#pragma once
#include "Eigen/Eigen"
#include "vector"

using namespace std;

struct TracesMesh {
    vector<unsigned int> traces_id; //vettore contenente gli id di ogni traccia
    vector<array<Eigen::Vector3d,2>> traces_vertices; //vettore contenente le coppie di vertici (definiti dalle rispettive coordinate) che delimitano la traccia
    vector<array<unsigned int,2>> traces_fracture; //vettore degli id delle 2 fratture che generano la traccia
    vector<double> traces_length; //vettore contenente le lunghezze di ogni traccia
    TracesMesh() = default;
    void addTrace(unsigned int& id, array<Eigen::Vector3d,2>& vertices, array<unsigned int,2>& fractures) {
        traces_id.push_back(id);
        traces_vertices.push_back(vertices);
        traces_fracture.push_back(fractures);
        traces_length.push_back((vertices[0] - vertices[1]).norm()); //calcolo lunghezza della traccia mediante la norma tra i due vertici che la delimitano
    }
};
