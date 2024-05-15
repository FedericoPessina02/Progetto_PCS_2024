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

    calculateSphere();
}

void Fracture::calculateSphere() {
    for (int i = 0; i <3; i++) {
        double average = 0;
        for (int j =  0; j < num_vertices; j++) {
            average += vertices(i, j);
        }
        average /= num_vertices;
        barycenter(i) = average;
    }

    double max_radius = 0;
    for (int i = 0; i < num_vertices; i++) {
        Eigen::Vector3d vertex = vertices.col(i) - barycenter;
        if (vertex.norm() > max_radius) {
            max_radius = vertex.norm();
        }
    }
    radius = max_radius;

}

void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {

}

PolygonalMesh generatePolygonalMesh() {
    PolygonalMesh output;
    return output;
}

