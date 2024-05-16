#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "PolygonalMesh.hpp"
#include "TracesMesh.hpp"
#include "Algorithms.hpp"

using namespace std;

Fracture::Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices) {
    id = _id;
    num_vertices = _num_vertices;
    vertices = _vertices;
    calculateNormalVector();
}

void Fracture::calculateNormalVector(){
    Eigen::Vector3d lato1;
    Eigen::Vector3d lato2;
    lato1 = vertices.col(1) - vertices.col(0);
    lato2 = vertices.col(2) - vertices.col(0);
    normal =lato1.cross(lato2);
    plane_d = normal.dot(vertices.col(0));
}

vector<Eigen::Vector3d> Fracture::calculateIntersectionsPoints(Eigen::Vector3d line, Eigen::Vector3d point) {
    vector<Eigen::Vector3d> result;
    for (unsigned int i = 0; i < num_vertices; i++) {
        Eigen::Vector3d lato;
        if (i == num_vertices-1) {
            lato = vertices.col(0) - vertices.col(i);
        } else {
            lato = vertices.col(i+1) - vertices.col(i);
        }
        if (lato.cross(line).norm() < 5*numeric_limits<double>::epsilon()) {
            continue; // il lato è parallelo alla traccia non ha senso cercare un'intersezione
        }
        Eigen::Matrix2d A;
        A << lato, -1*line;
        Eigen::Vector3d b = point - vertices.col(i);
        // riduco a due equazioni (significative)
        Eigen::Matrix2d Coeffs;
        Eigen::Vector2d b_values;
        if (abs(A.block(0, 0, 2, 2).determinant()) >= 5*numeric_limits<double>::epsilon()) {
            Coeffs = A.block(0, 0, 2, 2);
            b_values << b(0), b(1);
        } else if (abs(A.block(1, 0, 2, 2).determinant()) >= 5*numeric_limits<double>::epsilon()) {
            Coeffs = A.block(1, 0, 2, 2);
            b_values << b(1), b(2);
        } else {
            Coeffs << line(0), lato(0), line(2), lato(2); // controllare che funzioni
            b_values << b(0), b(2);
        }

        Eigen::Vector3d parameters = Coeffs.lu().solve(b_values);
        if (0<=parameters(0) && parameters(0)<=1) {
            Eigen::Vector3d point = vertices.col(i) + parameters(0)*lato;
            result.push_back(point);
        }
    }
    return result;
}


void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {
    //stiamo lavorando su una frattura, passiamo in input un'altra frattura per verificare
    //che ci sia intersezione e tramite referenza scrivo sulla mesh delle tracce
    array<Eigen::Vector3d,2> trace_verteces;
    unsigned int trace_id = mesh.traces_id.size();

    Eigen::Vector3d t = other.normal.cross(normal);
    Eigen::Matrix3d A;
    A.row(0) = normal;
    A.row(1) = other.normal;
    A.row(2) = t;
    if (abs(A.determinant()) < 5*numeric_limits<double>::epsilon()) {
        return;
    }
    Eigen::Vector3d b(plane_d, other.plane_d, 0);
    Eigen::Vector3d p = A.lu().solve(b); //gestire errore fattorizzazione

    vector<Eigen::Vector3d> punti_1 = calculateIntersectionsPoints(t, p);
    vector<Eigen::Vector3d> punti_2 = other.calculateIntersectionsPoints(t, p);
    // escludere caso di traccia impropria con i punti estremali e poi prendere i due interni
}

PolygonalMesh generatePolygonalMesh() {
    PolygonalMesh output;
    return output;
}

