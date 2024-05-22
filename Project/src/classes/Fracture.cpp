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

PolygonalMesh Fracture::generatePolygonalMesh(TracesMesh& traces) {
    PolygonalMesh mesh;
    // riempie la mesh poligonale con i punti e i lati iniziali e crea il primo poligono
    mesh.FractureId = id;
    for (unsigned int i = 0; i < num_vertices; i++) {
        mesh.CoordinateCell0Ds.push_back(vertices.col(i));
        mesh.IdCell0Ds.push_back(i);
        if (i<num_vertices-1) {
            mesh.IdCell1Ds.push_back(mesh.IdCell1Ds.size());
            mesh.VerticesCell1Ds.push_back(array<unsigned int,2> {i, i+1});
        } else {
            mesh.IdCell1Ds.push_back(mesh.IdCell1Ds.size());
            mesh.VerticesCell1Ds.push_back(array<unsigned int,2> {i, 0});
        }
    }
    mesh.IdCell2Ds.push_back(0);
    mesh.activatedPolygons.push_back(0);
    mesh.VerticesCell2Ds.push_back(mesh.IdCell0Ds);
    mesh.EdgesCell2Ds.push_back(mesh.IdCell1Ds);

    // taglio per tracce passanti
    for (unsigned int& trace_id: passant_traces) {
        Eigen::Vector3d direction = traces.traces_vertices[trace_id][1] - traces.traces_vertices[trace_id][0];
        Eigen::Vector3d application_point = traces.traces_vertices[trace_id][0];
        for (unsigned int& polygonId: mesh.activatedPolygons) {
            vector<Eigen::Vector3d> intersection_points;
            vector<unsigned int> intersection_starters;
            // scorro i lati del poligono per cercare i punti di intersezione
            for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygonId].size(); i++) {
                Eigen::Vector3d a = mesh.CoordinateCell0Ds[mesh.VerticesCell2Ds[polygonId][i]];
                Eigen::Vector3d b;
                if (i < mesh.VerticesCell2Ds[polygonId].size() - 1) {
                   b = mesh.CoordinateCell0Ds[mesh.VerticesCell2Ds[polygonId][i+1]];
                } else {
                   b = mesh.CoordinateCell0Ds[mesh.VerticesCell2Ds[polygonId][0]];
                }
                Eigen::Vector3d edge_direction = b-a;
                Eigen::Vector3d edge_application = a;
                if (edge_direction.cross(direction).norm() < 5*numeric_limits<double>::epsilon()) {
                    continue; // il lato è parallelo al segmento non ha senso cercare un'intersezione
                }
                Eigen::MatrixXd A;
                A.resize(3,2);
                A.col(0) = edge_direction;
                A.col(1) = -1*direction;
                Eigen::Vector3d coef = application_point - edge_application;
                Eigen::Vector2d parameters = A.colPivHouseholderQr().solve(coef);
                if (-5*numeric_limits<double>::epsilon()<=parameters(0) && parameters(0)<1+5*numeric_limits<double>::epsilon()) {
                    Eigen::Vector3d point = this->vertices.col(i) + parameters(0)*edge_direction;
                    intersection_points.push_back(point);
                    intersection_starters.push_back(mesh.VerticesCell2Ds[polygonId][i]);
                }
            }
            if (intersection_points.size() != 2) {
                continue;
            }
            array<unsigned int, 2> segment;
            segment[0] = mesh.addPoint(intersection_points[0]);
            segment[1] = mesh.addPoint(intersection_points[1]);
            array<unsigned int, 2> line_starters = {intersection_starters[0], intersection_starters[1]};
            Algorithms::cutPolygonBySegment(*this, mesh, polygonId, segment, line_starters);
        }
    }

    // taglio per tracce non passanti
    for (unsigned int& trace_id: internal_traces) {
        array<Eigen::Vector3d,2>& trace_vertices = traces.traces_vertices[trace_id];
        Eigen::Vector3d direction = traces.traces_vertices[trace_id][1] - traces.traces_vertices[trace_id][0];
        Eigen::Vector3d application_point = traces.traces_vertices[trace_id][0];
    }

    return mesh;
}
