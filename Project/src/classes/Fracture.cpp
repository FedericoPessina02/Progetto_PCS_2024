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
}

void Fracture::normalVector(){
    array<double,3> lato1;
    array<double,3> lato2;
    for (unsigned int i = 0; i < 3; i++)
    {
        lato1[i] = vertices(i,1) - vertices(i,0);
        lato2[i] = vertices(i,2) - vertices(i,0);
    }
    normal = crossProduct(lato1,lato2);
}


void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) { //stiamo lavorando su una frattura, passiamo in input un'altra frattura per verifi
                                                                  //care che ci sia intersezione e tramite referenza scrivo sulla mesh delle tracce
    array<double,3> t = crossProduct(other.normal,normal);

}

PolygonalMesh generatePolygonalMesh() {
    PolygonalMesh output;
    return output;
}

