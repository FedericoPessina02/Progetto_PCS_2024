#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"

using namespace std;

namespace Algorithms {

vector<Eigen::Vector3d> calculateIntersectionsPoints(Fracture& fracture, Eigen::Vector3d line, Eigen::Vector3d point) {
    vector<Eigen::Vector3d> result;
    for (unsigned int i = 0; i < fracture.num_vertices; i++) {
        Eigen::Vector3d lato;
        if (i == fracture.num_vertices-1) {
            lato = fracture.vertices.col(0) - fracture.vertices.col(i);
        } else {
            lato = fracture.vertices.col(i+1) - fracture.vertices.col(i);
        }
        if (lato.cross(line).norm() < 5*numeric_limits<double>::epsilon()) {
            continue; // il lato è parallelo alla traccia non ha senso cercare un'intersezione
        }
        Eigen::MatrixXd A;
        A.resize(3,2);
        A.col(0) = lato;
        A.col(1) = -1*line;
        Eigen::Vector3d b = point - fracture.vertices.col(i);
        Eigen::Vector2d parameters = A.colPivHouseholderQr().solve(b);
        if (0<=parameters(0) && parameters(0)<1) {
            Eigen::Vector3d point = fracture.vertices.col(i) + parameters(0)*lato;
            result.push_back(point);
        }
    }
    return result;
}

map<int, vector<Fracture>> assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number) {
    map<int, vector<Fracture>> id_to_fractures;
    // domain_borders = {x_min, x_max, y_min, y_max, z_min, z_max}
    // partitions_number è il numero di partizioni per dimensione (2-->8 partizioni, 3-->27 partizioni, ecc ecc)
    array<double, 3> chunk_size = {
        abs(domain_borders[1] - domain_borders[0]) / partitions_number,
        abs(domain_borders[3] - domain_borders[2]) / partitions_number,
        abs(domain_borders[5] - domain_borders[4]) / partitions_number
    };

    for(Fracture& el: fractures) {
        int x_partition_0 = floor(abs(el.vertices(0, 0) - domain_borders[0]) / chunk_size[0]);
        int y_partition_0 = floor(abs(el.vertices(1, 0) - domain_borders[2]) / chunk_size[1]);
        int z_partition_0 = floor(abs(el.vertices(2, 0) - domain_borders[4]) / chunk_size[2]);
        // converte l'id della partizione da base n=partitions_number a base 10 (id sequenziali) e aggiunge 1 (0 riservato alle fratture che sforano)
        int partition_id_vertex_0 = 1 + x_partition_0 + y_partition_0*partitions_number + z_partition_0*partitions_number*partitions_number;
        el.partition_id = partition_id_vertex_0;
        for(unsigned int i = 1; i < el.num_vertices; i++) {
            int x_partition = floor(abs(el.vertices(0, i) - domain_borders[0]) / chunk_size[0]);
            int y_partition = floor(abs(el.vertices(1, i) - domain_borders[2]) / chunk_size[1]);
            int z_partition = floor(abs(el.vertices(2, i) - domain_borders[4]) / chunk_size[2]);
            int partition_id_vertex_i = 1 + x_partition + y_partition*partitions_number + z_partition*partitions_number*partitions_number;
            if (partition_id_vertex_i != partition_id_vertex_0) {
                el.partition_id = 0;
                break;
            }
        }
        id_to_fractures[0].push_back(el);
    }
    fractures.clear();
    return id_to_fractures;
}

void cutTracesInsidePartition(vector<Fracture>& fractures, TracesMesh& mesh) {
    for (unsigned int i = 0; i < fractures.size() - 1; i ++) {
        for (unsigned int j = i+1; j<fractures.size(); j++) {
            Fracture& a = fractures[i];
            Fracture& b = fractures[j];
            if ((a.barycenter-b.barycenter).squaredNorm() < 2 * (a.radius + b.radius)) {
                a.generateTrace(b, mesh);
            }
        }
    }
}

void cutTracesOverlapping(vector<Fracture>& overlapping, vector<Fracture>& other_fractures, TracesMesh& mesh) {
    for (unsigned int i = 0; i < overlapping.size(); i ++) {
        for (unsigned int j = 0; j<other_fractures.size(); j++) {
            Fracture& a = overlapping[i];
            Fracture& b = other_fractures[j];
            if ((a.barycenter-b.barycenter).squaredNorm() < 2 * (a.radius + b.radius)) {
                a.generateTrace(b, mesh);
            }
        }
    }
}

void cutTraces(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh) {
    for (unsigned int id = 0; id < id_to_fractures.size(); id++) {
        if (id_to_fractures[id].size() == 0) {
            continue;
        }
        cutTracesInsidePartition(id_to_fractures[id], mesh);
        if (id != 0) {
            cutTracesOverlapping(id_to_fractures[0], id_to_fractures[id], mesh);
        }
    }
}

vector<PolygonalMesh> cutPolygonalMesh(map<int, vector<Fracture>>& id_to_fractures) {
    vector<PolygonalMesh> output;
    for (unsigned int id = 0; id < id_to_fractures.size(); id++) {
        for (Fracture& el: id_to_fractures[id]) {
            output.push_back(el.generatePolygonalMesh());
        }
    }
    return output;
}

void cutPolygonBySegment(PolygonalMesh& mesh, unsigned int polygonId, array<Eigen::Vector3d,2> segment) {
    // taglia il poligono (convesso) in due seguendo il segmento (passante per due suoi lati)
}

}
