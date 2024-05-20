#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "PolygonalMesh.hpp"
#include "TracesMesh.hpp"
#include "Algorithms.hpp"
#include "Utils.hpp"

using namespace std;

Fracture::Fracture(unsigned int& _id,unsigned int& _num_vertices, Eigen::MatrixXd& _vertices) {
    id = _id;
    num_vertices = _num_vertices;
    vertices = _vertices;
    calculateNormalVector();
    calculateSphere();
}

void Fracture::calculateSphere() {
    for (int i = 0; i <3; i++) {
        double average = 0;
        for (unsigned int j =  0; j < num_vertices; j++) {
            average += vertices(i, j);
        }
        average /= num_vertices;
        barycenter(i) = average;
    }

    double max_radius = 0;
    for (unsigned int i = 0; i < num_vertices; i++) {
        Eigen::Vector3d vertex = vertices.col(i) - barycenter;
        if (vertex.squaredNorm() > max_radius) {
            max_radius = vertex.squaredNorm();
        }
    }
    radius = max_radius;
}

void Fracture::calculateNormalVector(){
    Eigen::Vector3d lato1;
    Eigen::Vector3d lato2;
    lato1 = vertices.col(1) - vertices.col(0);
    lato2 = vertices.col(2) - vertices.col(0);
    normal =lato1.cross(lato2);
    plane_d = normal.dot(vertices.col(0));
}

void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {
    //stiamo lavorando su una frattura, passiamo in input un'altra frattura per verificare
    //che ci sia intersezione e tramite referenza scrivo sulla mesh delle tracce
    array<Eigen::Vector3d,2> trace_verteces;
    unsigned int trace_id = mesh.traces_id.size();

    Eigen::Vector3d t = normal.cross(other.normal);
    Eigen::Matrix3d A;
    A.row(0) = normal;
    A.row(1) = other.normal;
    A.row(2) = t;
    if (abs(A.determinant()) < 5*numeric_limits<double>::epsilon()) {
        return;
    }
    Eigen::Vector3d b(plane_d, other.plane_d, 0);
    Eigen::Vector3d p = A.lu().solve(b); //gestire errore fattorizzazione

    vector<Eigen::Vector3d> punti_1 = Algorithms::calculateIntersectionsPoints(*this, t, p);
    vector<Eigen::Vector3d> punti_2 = Algorithms::calculateIntersectionsPoints(other, t, p);

    if (punti_1.size()+punti_2.size()!=4) {
        if (punti_1.size()+punti_2.size()>4) {
            cerr << "More than 4 points of intersections";
            return;
        }
        if (punti_1.size() == 1 || punti_2.size() == 1) {
            cerr << "1 punto di intersezione" << endl;
            return;
        }
        return;
    }

    vector<Eigen::Vector3d> punti_distinti = Utils::calculateDistinctPoints(punti_1, punti_2);

    if (punti_distinti.size() == 2) {
        array<Eigen::Vector3d, 2> punti_distinti_array;
        copy_n(make_move_iterator(punti_distinti.begin()), 2, punti_distinti_array.begin());
        array<unsigned int, 2> fractures_id = {id, other.id};
        mesh.addTrace(trace_id, punti_distinti_array, fractures_id);
        passant_traces.push_back(id);
        other.passant_traces.push_back(id);
        return;
    }
    if (punti_distinti.size() == 3) {
        Eigen::Vector3d doppione;
        for (Eigen::Vector3d& i: punti_distinti) {
            unsigned int counter = 0;
            for (Eigen::Vector3d& a: punti_1) {
                if ((i-a).squaredNorm() < 8*numeric_limits<double>::epsilon()) {
                    counter += 1;
                }
            }
            for (Eigen::Vector3d& a: punti_2) {
                if ((i-a).squaredNorm() < 8*numeric_limits<double>::epsilon()) {
                    counter += 1;
                }
            }
            if (counter == 2) {
                doppione = i;
                break;
            }
        }

        double min_length = numeric_limits<double>::max();
        array<Eigen::Vector3d, 2> punti_distinti_array;
        for (Eigen::Vector3d& i: punti_distinti) {
            if (i == doppione) {
                continue;
            }
            double new_length = (doppione-i).squaredNorm();
            if (new_length < min_length) {
                min_length = new_length;
                punti_distinti_array = {doppione, i};
            }
        }
        array<unsigned int, 2> fractures_id = {id, other.id};
        mesh.addTrace(trace_id, punti_distinti_array, fractures_id);

        double punti_1_length = (punti_1[0] - punti_1[1]).squaredNorm();
        if (abs(min_length-punti_1_length) < 8*numeric_limits<double>::epsilon()) {
            passant_traces.push_back(trace_id);
            other.internal_traces.push_back(trace_id);
        } else {
            internal_traces.push_back(trace_id);
            other.passant_traces.push_back(trace_id);
        }
        return;
    }

    vector<Eigen::Vector3d>& punti = punti_1;
    double distanza = 0;
    Eigen::Vector3d estremo;
    for (unsigned int i = 0; i < punti_1.size(); i++) {
        for (unsigned int j = 0; j < punti_2.size(); j++) {
            if ((punti_1[i]-punti_2[j]).squaredNorm() > distanza) {
                estremo = punti_1[i];
            }
        }
    }
    if ((punti_2[0]-punti_2[1]).squaredNorm() > distanza) {
        estremo = punti_2[0];
        punti = punti_2;
    }

    Eigen::Vector3d punto_vicini;
    distanza = numeric_limits<double>::max();
    for (Eigen::Vector3d& i: punti_distinti) {
        if (i == estremo) {
            continue;
        }
        if ((i-estremo).squaredNorm() < distanza) {
            punto_vicini = i;
            distanza = (i-estremo).squaredNorm();
        }
    }

    if (find(punti.begin(), punti.end(), punto_vicini) != punti.end()) {
        return;
    }

    Eigen::Vector3d punto_meno_vicini;
    distanza = numeric_limits<double>::max();
    for (Eigen::Vector3d& i: punti_distinti) {
        if (i == estremo || i == punto_vicini) {
            continue;
        }
        if ((i-estremo).squaredNorm() < distanza) {
            punto_meno_vicini = i;
            distanza = (i-estremo).squaredNorm();
        }
    }

    vector<Eigen::Vector3d> punti_distinti_vector = {punto_vicini, punto_meno_vicini};
    array<Eigen::Vector3d, 2> punti_distinti_array = {punto_vicini, punto_meno_vicini};
    array<unsigned int, 2> fractures_id = {id, other.id};
    mesh.addTrace(trace_id, punti_distinti_array, fractures_id);

    if (Utils::compareSegments(punti_distinti_vector, punti_1)) {
        passant_traces.push_back(trace_id);
    } else {
        internal_traces.push_back(trace_id);
    }

    if (Utils::compareSegments(punti_distinti_vector, punti_2)) {
        other.passant_traces.push_back(trace_id);
    } else {
        other.internal_traces.push_back(trace_id);
    }
}

PolygonalMesh Fracture::generatePolygonalMesh() {
    PolygonalMesh mesh;
    // riempie la mesh poligonale con i punti iniziali e crea il primo poligono
    mesh.FractureId = id;
    for (unsigned int i = 0; i < num_vertices; i++) {
        mesh.CoordinateCell0Ds.push_back(vertices.col(i));
        mesh.IdCell0Ds.push_back(i+1);
    }
    mesh.IdCell2Ds.push_back(1);
    mesh.VerticesCell2Ds.push_back(mesh.IdCell0Ds);
    return mesh;
}