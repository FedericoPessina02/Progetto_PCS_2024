#include "Fracture.hpp"
#include "Eigen/Eigen"
#include "vector"
#include "PolygonalMesh.hpp"
#include "TracesMesh.hpp"

using namespace std;

Fracture::Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices) {
    id = _id;
    num_vertices = _num_vertices;
    vertices = _vertices;
}

void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {

}

PolygonalMesh generatePolygonalMesh() {
    PolygonalMesh output;
    return output;
}

