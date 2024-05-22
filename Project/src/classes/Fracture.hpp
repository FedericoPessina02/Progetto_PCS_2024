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
    unsigned int partition_id = 0; //partizione dello spazio globale in cui è contenuta la frattura
    Eigen::MatrixXd vertices; //matrice contenente le coordinate dei vertici della frattura: colonne= vertici ; righe= x,y,z
    vector<unsigned int> internal_traces;
    vector<unsigned int> passant_traces;
    Eigen::Vector3d normal; //normale del piano contenente la frattura
    double plane_d;
    Eigen::Vector3d barycenter; //centroide della frattura
    double radius; //raggio della sfera contenente la frattura (NB! è al quadrato per evitare gli alti costi computazionali delle radici)
    Fracture() = default;
    Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices);
    void generateTrace(Fracture& other, TracesMesh& mesh);
    void calculateNormalVector();
    void calculateSphere();
    vector<Eigen::Vector3d> calculateIntersectionsPoints(Eigen::Vector3d line, Eigen::Vector3d point);
    PolygonalMesh generatePolygonalMesh();
};

