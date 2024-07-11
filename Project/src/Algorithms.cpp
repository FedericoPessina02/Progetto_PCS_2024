#include "Fracture.hpp"
#include <iostream>
#include <fstream>
#include <thread>
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
        if (lato.cross(line).norm() < Utils::tol_coeff*numeric_limits<double>::epsilon()) {
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
        if (-Utils::tol_coeff*numeric_limits<double>::epsilon()<parameters(0) && parameters(0)<1-Utils::tol_coeff*numeric_limits<double>::epsilon()) {
            // verifico che l'intersezione con il lato sia una combinazione convessa, ossia che effettivamente il punto di intersezione
            // cada all'interno del lato
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
        for(unsigned int i = 1; i < el.num_vertices; i++) {
            //verifica che tutti i vertici della frattura siano nella stessa partizione
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
        // la frattura viene aggiunta alla mappa nel vettore associato al suo id
        id_to_fractures[el.partition_id].push_back(el);
    }
    // il vettore di fratture iniziali viene cancellato dalla memoria per liberare spazio
    fractures.clear(); // <- rimuovo gli elementi dal vettore
    fractures.shrink_to_fit(); // <- libero la memoria riservata
    return id_to_fractures;
}

void cutTracesInsidePartition(vector<Fracture>& fractures, TracesMesh& mesh) {
    if (fractures.size() == 0) {
        return;
    }
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
    // funzione di smistamento per gestire il taglio delle tracce
    // analizza tutte le partizioni
    for (unsigned int id = 0; id < id_to_fractures.size(); id++) {
        if (id_to_fractures[id].size() == 0) {
            continue;
        }
        cutTracesInsidePartition(id_to_fractures[id], mesh); // calcolo delle tracce generate all'interno della partizione
        if (id != 0) {
            // se sto analizzando le fratture di una partizione ben definita, devo anche confrontarle con le fratture con label 0
            cutTracesOverlapping(id_to_fractures[0], id_to_fractures[id], mesh);
        }
    }
}

// void cutTracesMultithread(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh) {
//     auto standardCut = [](vector<Fracture>& fractures, TracesMesh& mesh) {
//         cutTracesInsidePartition(fractures, mesh);
//     };
//     auto overlapCut = [](vector<Fracture>& fractures_a, vector<Fracture>& fractures_b, TracesMesh& mesh) {
//         cutTracesOverlapping(fractures_a, fractures_b, mesh);
//     };

//     vector<thread> processes;

//     for (int id_zone = 0; id_zone < id_to_fractures.size() / 4; id_zone++) {
//         processes.clear();
//         for (int i = 0; i < 4; i++) {
//             processes.push_back(thread(standardCut, ref(id_to_fractures[id_zone+i]), ref(mesh)));
//         }
//         for(auto& th : processes){
//             th.join();
//         }

//         processes.clear();
//         for (int i = 0; i < 4; i++) {
//             if (id_zone + 1 == 0) {
//                 continue; // il caso 0 contro 0 è gia un caso standard trattato in precedenza
//             }
//             processes.push_back(thread(overlapCut, ref(id_to_fractures[0]), ref(id_to_fractures[id_zone+i]), ref(mesh)));
//         }
//         for(auto& th : processes){
//             th.join();
//         }
//     }

//     int position = id_to_fractures.size() / 4;
//     int remain = id_to_fractures.size() % 4;
//     processes.clear();
//     for (int i = 0; i < remain; i++) {
//         processes.push_back(thread(standardCut, ref(id_to_fractures[position+i]), ref(mesh)));
//         if (position+i == 0) {
//             continue;
//         }
//         processes.push_back(thread(overlapCut, ref(id_to_fractures[0]), ref(id_to_fractures[position+i]), ref(mesh)));
//     }

//     for(auto& th : processes){
//         th.join();
//     }
// }

vector<PolygonalMesh> cutPolygonalMesh(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& traces_mesh) {
    // funzione di smistamente per avviare il taglio delle mesh poligonali associate a ogni frattura
    // successivamente le memorizza nel vettore di output
    vector<PolygonalMesh> output;
    for (unsigned int fracture_id = 0; fracture_id < id_to_fractures.size(); fracture_id++) {
        for (Fracture& el: id_to_fractures[fracture_id]) {
            output.push_back(el.generatePolygonalMesh(traces_mesh));
        }
    }
    return output;
}

vector<PolygonalMesh> cutPolygonalMeshMultithread(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& traces_mesh) {
    // lambda function da chiamare nei thread
    auto createMesh = [](vector<PolygonalMesh>& output, Fracture& fracture, TracesMesh& mesh) {
        output.push_back(fracture.generatePolygonalMesh(mesh));
    };

    // vettore contenente le mesh generate
    vector<PolygonalMesh> output;
    // vettore contenente i thread aperti (sono al massimo 4)
    vector<thread> processes;

    // le mesh vengono calcolate a gruppi di 4, uno per thread
    // il risultato è analogo a cutPolygonalMesh, ma ne vengono eseguiti 4 in parallelo (i calcoli sono distinti non ci sono problemi di sincronizzazione)
    for (unsigned int fracture_id = 0; fracture_id < id_to_fractures.size(); fracture_id++) {
        for (int id_zone = 0; id_zone < id_to_fractures[fracture_id].size() / 4; id_zone++) {
            // svuoto il vettore contenente i thread
            processes.clear();
            // avvio i 4 thread assegnando ad ognuno una frattura da tagliare
            for (int i = 0; i < 4; i++) {
                processes.push_back(thread(createMesh, ref(output), ref(id_to_fractures[fracture_id][4*id_zone+i]), ref(traces_mesh)));
            }
            // aspetto che tutti e 4 abbiano finito prima di ripetere il processo
            for(auto& th : processes){
                th.join();
            }
        }

        // il codice è analogo a prima ma si occupa delle ultime fratture (se il totale non è un multiplo di 4 ne processa 1/2/3)
        int remain = id_to_fractures[fracture_id].size() % 4;
        int position = id_to_fractures[fracture_id].size() - remain;
        processes.clear();
        for (int i = 0; i < remain; i++) {
            processes.push_back(thread(createMesh, ref(output), ref(id_to_fractures[fracture_id][position+i]), ref(traces_mesh)));
        }

        for(auto& th : processes){
            th.join();
        }
    }
    return output;
}

void cutPolygonBySegment(Fracture& fracture, PolygonalMesh& mesh, unsigned int polygonId, vector<unsigned int> total_points, array<unsigned int,2> segment) {
    // taglia il poligono in due sottopoligoni seguendo il segmento avente come estremi due vertici sul perimetro del poligono
    // come prima cosa creo i due vettori che conterranno gli id dei vertici dei due sottopoligoni risultanti
    vector<unsigned int> polygon_a_vertices;
    vector<unsigned int> polygon_b_vertices;
    unsigned int CUT_ALG = 2; // stabilisce quale algoritmo usare
    if (CUT_ALG == 1) {
        /* ALGORITMO 1
        Questo algoritmo sfrutta i prodotti vettoriali per stabilire a quale sottopoligono appartengono i vari vertici
        Scorrendo uno ad uno tutti i vertici faccio il prodotto vettoriale tra il segmento di taglio e la congiungente del vertice in questione con un punto del segmento
        il vettore ottenuto mi permette di discriminare il sottopoligono di appartenenza */
        Eigen::Vector3d cut_line = mesh.CoordinateCell0Ds[segment[1]] - mesh.CoordinateCell0Ds[segment[0]];
        for (unsigned int& vertex_id: total_points) {
            if (vertex_id == segment[0] || vertex_id == segment[1]) {
                // se il punto in questione appartiene al segmento di taglio esso sarà necessariamente in entrambi i sottopoligoni risultanti
                polygon_a_vertices.push_back(vertex_id);
                polygon_b_vertices.push_back(vertex_id);
                continue;
            }
            Eigen::Vector3d link_line = mesh.CoordinateCell0Ds[vertex_id] - mesh.CoordinateCell0Ds[segment[0]];
            Eigen::Vector3d product_line = cut_line.cross(link_line);
            double evaluation_coef = fracture.normal.dot(product_line) + fracture.plane_d;
            // se sotituisco il prodotto vettoriale nell'equazione planare posso discriminare a quale sottopoligono appartiene il punto guardando il segno
            // la scelta di associare al primo sottopoligono i valori positivi è totalmente arbitraria...
            if (evaluation_coef > Utils::tol_coeff*numeric_limits<double>::epsilon()) {
                polygon_a_vertices.push_back(vertex_id);
            } else if(evaluation_coef < -Utils::tol_coeff*numeric_limits<double>::epsilon()) {
                polygon_b_vertices.push_back(vertex_id);
            }
        }
    } else if (CUT_ALG == 2) {
        /* ALGORITMO 2
        Questo algoritmo sfrutta la lettura sequenziale del vettore ordinato in senso antiorario dei punti
        parte con l'assegnare i punti al primo sottopoligono, e quando si imbatte in un punto appartenente al segmento di taglio
        inverte la variabile booleana e comincia a memorizzare i punti successivi nell'altro sottopoligono (i punti di taglio vengono oviamente memorizzati in entrambi i sottopoligoni */
        // l'algoritmo per il resto è molto semplice
        bool polygon_flag = true;
        for (unsigned int& vertex_id: total_points) {
            if (vertex_id == segment[0] || vertex_id == segment[1]) {
                polygon_a_vertices.push_back(vertex_id);
                polygon_b_vertices.push_back(vertex_id);
                polygon_flag = !polygon_flag;
            } else {
                if (polygon_flag) {
                    polygon_a_vertices.push_back(vertex_id);
                } else {
                    polygon_b_vertices.push_back(vertex_id);
                }
            }
        }
    }

    //ora che ho i due vettori contenenti i punti distribuiti tra i due sottopoligoni devo aggiornare la mesh
    // come prima cosa creo due nuovi poligoni e vi assegno le rispettive liste di punti
    // notare che gli algoritmi non modificano l'ordine dei punti --> l'integrità della mesh di partenza è preservata -sempre-
    unsigned int id_1 = mesh.IdCell2Ds.size();
    mesh.IdCell2Ds.push_back(id_1);
    mesh.VerticesCell2Ds.push_back(polygon_a_vertices);
    unsigned int id_2 = mesh.IdCell2Ds.size();
    mesh.IdCell2Ds.push_back(id_2);
    mesh.VerticesCell2Ds.push_back(polygon_b_vertices);

    // il poligono padre viene cancellato dall'elenco dei poligoni attivi e al contempo memorizzo i due poligoni nuovi come attivi
    mesh.activatedPolygons.erase(remove(mesh.activatedPolygons.begin(), mesh.activatedPolygons.end(), polygonId), mesh.activatedPolygons.end());
    mesh.activatedPolygons.push_back(id_1);
    mesh.activatedPolygons.push_back(id_2);

    // aggiorno ora i lati della mesh
    // aggiungo il primo vertice di nuovo in fondo per poter fare la lettura ciclica dei lati senza dover distinguere casi particolari
    //     posso modificarlo perché tanto ho copiato il contenuto in mesh.VerticesCell2Ds, quindi sono effettivamente due oggetti separati
    polygon_a_vertices.push_back(polygon_a_vertices[0]);
    vector<unsigned int> polygon_a_edges;
    // scorro i vertici del poligono e salvo il lato (se esiste già mi restituisce l'id del lato preesistente)
    for (unsigned int i = 0; i < polygon_a_vertices.size() - 1; i++) {
        unsigned int edge_id = mesh.addEdge(polygon_a_vertices[i], polygon_a_vertices[i+1]);
        polygon_a_edges.push_back(edge_id);
    }
    mesh.EdgesCell2Ds.push_back(polygon_a_edges);

    // faccio lo stesso procedimento per il secondo sottopoligono
    polygon_b_vertices.push_back(polygon_b_vertices[0]);
    vector<unsigned int> polygon_b_edges;
    for (unsigned int i = 0; i < polygon_b_vertices.size() - 1; i++) {
        unsigned int edge_id = mesh.addEdge(polygon_b_vertices[i], polygon_b_vertices[i+1]);
        polygon_b_edges.push_back(edge_id);
    }
    mesh.EdgesCell2Ds.push_back(polygon_b_edges);
}

void ordinaTracce(map<int, vector<Fracture>>& id_to_fractures, TracesMesh& mesh, string nome_file){
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
                // BubbleSort sembra vincere? Alla fine gli insiemi da ordinare sono sempre relativamente piccoli...
                // sort(vec.begin(), vec.end(), compareByValueDescending);
                // sort(vect.begin(), vect.end(), compareByValueDescending);
                // Utils::MergeSort(vec);
                // Utils::MergeSort(vect);
                Utils::BubbleSort(vec);
                Utils::BubbleSort(vect);

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


