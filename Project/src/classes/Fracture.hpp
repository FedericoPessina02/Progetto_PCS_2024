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
    Eigen::Vector3d normal; // normale del piano
    double plane_d; // coefficiente dell'equazione del piano
    Eigen::Vector3d barycenter; // baricentro della frattura
    double radius; // raggio al quadrato della sfera che contiene la frattura
    Fracture() = default;
    Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices);
    void calculateNormalVector();
    void calculateSphere();
    void generateTrace(Fracture& other, TracesMesh& mesh);
    PolygonalMesh generatePolygonalMesh(TracesMesh& traces);
};

