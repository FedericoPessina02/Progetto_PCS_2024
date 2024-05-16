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
    vector<unsigned int> passant_traces;
    Fracture() = default;
    Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices);
    void generateTrace(Fracture& other, TracesMesh& mesh);
    void calculateNormalVector();
    vector<Eigen::Vector3d> calculateIntersectionsPoints(Eigen::Vector3d line, Eigen::Vector3d point);
    PolygonalMesh generatePolygonalMesh();
    Eigen::Vector3d normal; //mi salvo la normale del piano contenente il poligono
    double plane_d;
};

