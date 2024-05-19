#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"

using namespace std;

namespace Algorithms {

void assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number) {
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
    }
}

void cutTraces(vector<Fracture>& fractures, TracesMesh& mesh) {
    for (unsigned int i = 0; i < fractures.size() - 1; i ++) {
        for (unsigned int j = i+1; j<fractures.size(); j++) {
            Fracture& a = fractures[i];
            Fracture& b = fractures[j];
            if (!(a.partition_id == b.partition_id || a.partition_id == 0 || b.partition_id == 0)) {
                continue; // vuol dire che le due fratture appartengono a partizioni ben precise e diverse
            }
            if ((a.barycenter-b.barycenter).squaredNorm() < 2 * (a.radius + b.radius)) {
                a.generateTrace(b, mesh);
            }
        }
    }
}

}
