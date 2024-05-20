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
    unsigned int partition_id = 0;
    Eigen::MatrixXd vertices;
    vector<unsigned int> internal_traces;
    vector<unsigned int> passant_traces;
    Eigen::Vector3d normal; //mi salvo la normale del piano contenente il poligono
    double plane_d;
    Eigen::Vector3d barycenter;
    double radius;
    Fracture() = default;
    Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices);
    void calculateNormalVector();
    void calculateSphere();
    void generateTrace(Fracture& other, TracesMesh& mesh);
    PolygonalMesh generatePolygonalMesh();
};

