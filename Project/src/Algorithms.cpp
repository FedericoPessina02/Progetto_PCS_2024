#include "Fracture.hpp"
#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "vector"
#include "Utils.hpp"

using namespace std;

namespace Algorithms {

vector<Eigen::Vector3d> calculateIntersectionsPoints(Fracture& fracture, Eigen::Vector3d line, Eigen::Vector3d point) { 
  //calcolo del punto di intersezione tra i lati della frattura e la retta di direzione 'line' e application point 'point' passati in input
  vector<Eigen::Vector3d> result;
    for (unsigned int i = 0; i < fracture.num_vertices; i++) {
        Eigen::Vector3d lato;
        //calcolo componenti dei lati della frattura, mediante differenza tra le coordinate di due vertici successivi
        //(NB! nel caso dell' ultimo vertice, useremo la differenza tra il primo e l'ultimo)
        if (i == fracture.num_vertices-1) {
            lato = fracture.vertices.col(0) - fracture.vertices.col(i);
        } else {
            lato = fracture.vertices.col(i+1) - fracture.vertices.col(i);
        }
        //verifica che il lato non sia parallelo alla traccia (in caso contrario, non c'è intersezione)
        if (lato.cross(line).norm() < 5*numeric_limits<double>::epsilon()) {
            continue;
        }
        //risoluzione del sistema lineare mediante la fattorizzazione QR (utile per matrici generiche mxn, avente il minor costo
        //computazionale ovvero (2n^3)/3  )
        Eigen::MatrixXd A;
        A.resize(3,2);
        A.col(0) = lato;
        A.col(1) = -1*line;
        Eigen::Vector3d b = point - fracture.vertices.col(i);
        Eigen::Vector2d parameters = A.colPivHouseholderQr().solve(b);
        if (-10*numeric_limits<double>::epsilon()<parameters(0) && parameters(0)<1-10*numeric_limits<double>::epsilon()) {
            Eigen::Vector3d point = fracture.vertices.col(i) + parameters(0)*lato;
            result.push_back(point);
        }
    }
    return result;
}

map<int, vector<Fracture>> assignPartition(vector<Fracture>& fractures, array<double, 6>& domain_borders, const int partitions_number) {
    //assegna ogni frattura ad una partizione dello spazio totale. Nel caso una frattura abbia vertici in più partizioni, viene assegnato il valore 0

    map<int, vector<Fracture>> id_to_fractures;
    /* domain_borders = {x_min, x_max, y_min, y_max, z_min, z_max}
       partitions_number è il numero di partizioni per dimensione (2-->8 partizioni, 3-->27 partizioni, ecc ecc)*/
    array<double, 3> chunk_size = {//calcola per x,y,z la dimensione del lato delle partizioni
        abs(domain_borders[1] - domain_borders[0]) / partitions_number,
        abs(domain_borders[3] - domain_borders[2]) / partitions_number,
        abs(domain_borders[5] - domain_borders[4]) / partitions_number
    };

    for(Fracture& el: fractures) {//dice in quale partizione dello spazio appartiene ciascuna frattura
        int x_partition_0 = floor(abs(el.vertices(0, 0) - domain_borders[0]) / chunk_size[0]);
        int y_partition_0 = floor(abs(el.vertices(1, 0) - domain_borders[2]) / chunk_size[1]);
        int z_partition_0 = floor(abs(el.vertices(2, 0) - domain_borders[4]) / chunk_size[2]);
        // converte l'id della partizione da base n=partitions_number a base 10 (id sequenziali) e aggiunge 1 (0 riservato alle fratture che sforano)
        int partition_id_vertex_0 = 1 + x_partition_0 + y_partition_0*partitions_number + z_partition_0*partitions_number*partitions_number;
        el.partition_id = partition_id_vertex_0;
        for(unsigned int i = 1; i < el.num_vertices; i++) {//verifica che tutti i vertici della frattura siano nella stessa partizione
            //altrimenti viene assegnato il valore 0
            int x_partition = floor(abs(el.vertices(0, i) - domain_borders[0]) / chunk_size[0]);
            int y_partition = floor(abs(el.vertices(1, i) - domain_borders[2]) / chunk_size[1]);
            int z_partition = floor(abs(el.vertices(2, i) - domain_borders[4]) / chunk_size[2]);
            int partition_id_vertex_i = 1 + x_partition + y_partition*partitions_number + z_partition*partitions_number*partitions_number;
            if (partition_id_vertex_i != partition_id_vertex_0) {
                el.partition_id = 0;
                break;
            }
        }
        id_to_fractures[el.partition_id].push_back(el);
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

vector<PolygonalMesh> cutPolygonalMesh(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& traces_mesh) {
    vector<PolygonalMesh> output;
    for (unsigned int id = 0; id < id_to_fractures.size(); id++) {
        for (Fracture& el: id_to_fractures[id]) {
            output.push_back(el.generatePolygonalMesh(traces_mesh));
        }
    }
    return output;
}

void cutPolygonBySegment(Fracture& fracture, PolygonalMesh& mesh, unsigned int polygonId, array<unsigned int,2> segment, array<unsigned int,2> intersection_starters) {
    // taglia il poligono (convesso) in due seguendo il segmento (passante per due suoi lati)
    vector<unsigned int> polygon_a_vertices;
    vector<unsigned int> polygon_b_vertices;
    Eigen::Vector3d cut_line = mesh.CoordinateCell0Ds[segment[1]] - mesh.CoordinateCell0Ds[segment[0]];
    for (unsigned int& vertex_id: mesh.VerticesCell2Ds[polygonId]) {
        Eigen::Vector3d link_line = mesh.CoordinateCell0Ds[vertex_id] - mesh.CoordinateCell0Ds[segment[0]];
        Eigen::Vector3d product_line = cut_line.cross(link_line);
        double evaluation_coef = fracture.normal.dot(product_line) + fracture.plane_d;
        if (evaluation_coef > 0) {
            polygon_a_vertices.push_back(vertex_id);
        } else if(evaluation_coef < 0) {
            polygon_b_vertices.push_back(vertex_id);
        } else {
            cerr << "Cut segment is parallel to an edge of the polygon";
            return;
        }
        if (vertex_id == intersection_starters[0]) {
            polygon_a_vertices.push_back(segment[0]);
            polygon_b_vertices.push_back(segment[0]);
        }
        if (vertex_id == intersection_starters[1]) {
            polygon_a_vertices.push_back(segment[1]);
            polygon_b_vertices.push_back(segment[1]);
        }
    }
    // creo un nuovo poligono assegnandogli gli estremi dell'altro poligono risultante
    unsigned int id_1 = mesh.IdCell2Ds.size();
    mesh.IdCell2Ds.push_back(id_1);
    mesh.VerticesCell2Ds.push_back(polygon_a_vertices);

    unsigned int id_2 = mesh.IdCell2Ds.size();
    mesh.IdCell2Ds.push_back(id_2);
    mesh.VerticesCell2Ds.push_back(polygon_b_vertices);

    mesh.activatedPolygons.erase(remove(mesh.activatedPolygons.begin(), mesh.activatedPolygons.end(), polygonId), mesh.activatedPolygons.end());
    mesh.activatedPolygons.push_back(id_1);
    mesh.activatedPolygons.push_back(id_2);

    polygon_a_vertices.push_back(polygon_a_vertices[0]);
    vector<unsigned int> polygon_a_edges;
    for (unsigned int i = 0; i < polygon_a_vertices.size() - 1; i++) {
        unsigned int edge_id = mesh.addEdge(polygon_a_vertices[i], polygon_a_vertices[i+1]);
        polygon_a_edges.push_back(edge_id);
    }
    mesh.EdgesCell2Ds.push_back(polygon_a_edges);

    polygon_b_vertices.push_back(polygon_b_vertices[0]);
    vector<unsigned int> polygon_b_edges;
    for (unsigned int i = 0; i < polygon_b_vertices.size() - 1; i++) {
        unsigned int edge_id = mesh.addEdge(polygon_b_vertices[i], polygon_b_vertices[i+1]);
        polygon_b_edges.push_back(edge_id);
    }
    mesh.EdgesCell2Ds.push_back(polygon_b_edges);
}

void ordinaFract(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh, string nome_file){
    ofstream ofs(nome_file);
    if (! ofs.is_open()){
        cerr<< "errore di apertura del file di output \n";
    }

    // Funzione che compara, utilizzata per ordinare in ordine decrescente
    auto compareByValueDescending = [](const pair<int, double>& a, const pair<int, double>& b) {
        return a.second > b.second;
    };
    for (unsigned int partition_id = 0; partition_id < id_to_fractures.size(); partition_id++) {
        if (id_to_fractures[partition_id].size() != 0) {
            for (Fracture& fract : id_to_fractures[partition_id]){
                ofs<<"# FractureId; NumTraces"<<'\n';
                ofs<<fract.id<<" ; ";
                ofs<<fract.internal_traces.size()+fract.passant_traces.size()<<'\n';
                ofs<<"# TraceId; Tips; Length"<<'\n';

                map<int,double> lunghezze_passanti; //mappa con chiave=id_traccia, valore=lunghezza_traccia
                map<int,double> lunghezze_interne;

                for (unsigned int& elem :fract.internal_traces){
                    lunghezze_interne[elem]=mesh.traces_length[elem];
                }

                for (unsigned int& elem :fract.passant_traces){ //Id di ogni traccia passante
                    lunghezze_passanti[elem]=mesh.traces_length[elem];
                }

                //crea un vettore di coppie (vec) contenente tutte le coppie chiave-valore presenti nella mappa lunghezze_passanti.
                vector<pair<int, double>> vec(lunghezze_passanti.begin(), lunghezze_passanti.end());
                //crea un vettore di coppie (vect) contenente tutte le coppie chiave-valore presenti nella mappa lunghezze_interne.
                vector<pair<int, double>> vect(lunghezze_interne.begin(), lunghezze_interne.end());

                // Ordinamento del vettore di coppie in base ai valori
                /* usare bubblesort o mergesort */
                sort(vec.begin(), vec.end(), compareByValueDescending);
                sort(vect.begin(), vect.end(), compareByValueDescending);

                for (const auto& pair : vect) {
                    ofs<< pair.first <<" ; ";
                    ofs<<"true"<<" ; ";
                    ofs<< pair.second<<'\n';
                }
                for (const auto& pair : vec) {
                    ofs<< pair.first <<" ; ";
                    ofs<<"false"<<" ; ";
                    ofs<< pair.second<<'\n';
                }
            }
        }
    }
    ofs.close();
}


}


