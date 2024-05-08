#pragma once
#include "Eigen/Eigen"
#include "vector"
#include "PolygonalMesh.hpp"
#include "TracesMesh.hpp"

using namespace std;

class Fracture {
public:
    unsigned int id;
    unsigned int num_vertices;
    Eigen::MatrixXd vertices;
    vector<unsigned int> internal_traces;
    Fracture() = default;
    Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices);
    void generateTrace(Fracture& other, TracesMesh& mesh);
    PolygonalMesh generatePolygonalMesh();
private:
    double barycenter;
    double radius;
};

