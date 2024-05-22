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
    //calcolo coordinate del centroide, mediante la media delle coordinate dei vertici rispetto a ciascuna coordinata
    for (int i = 0; i <3; i++) {
        double average = 0;
        for (unsigned int j =  0; j < num_vertices; j++) {
            average += vertices(i, j);
        }
        average /= num_vertices;
        barycenter(i) = average;
    }
    
    //calcolo della massima distanza tra il centroide e i vertici della frattura
    //NB! La distanza è al quadrato per evitare di calcolare la radice, che computazionalmente ha costo elevato
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
    //calcolo della normale al piano contenente la frattura. (Il piano è quello passante per i primi 3 vertici del poligono)
    Eigen::Vector3d lato1;
    Eigen::Vector3d lato2;
    lato1 = vertices.col(1) - vertices.col(0);
    lato2 = vertices.col(2) - vertices.col(0);
    normal =lato1.cross(lato2);
    plane_d = normal.dot(vertices.col(0));
}

vector<Eigen::Vector3d> Fracture::calculateIntersectionsPoints(Eigen::Vector3d line, Eigen::Vector3d point) {
    //calcolo del punto di intersezione tra i lati della frattura e la retta di direzione 'line' e application point 'point' passati in input

    vector<Eigen::Vector3d> result; //vettore contenente tutti i punti di intersezione tra la retta e i lati
    for (unsigned int i = 0; i < num_vertices; i++) {
        Eigen::Vector3d lato;

        //calcolo componenti dei lati della frattura, mediante differenza tra le coordinate di due vertici successivi
        //(NB! nel caso dell' ultimo vertice, useremo la differenza tra il primo e l'ultimo)
        if (i == num_vertices-1) {
            lato = vertices.col(0) - vertices.col(i);
        } else {
            lato = vertices.col(i+1) - vertices.col(i);
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
        Eigen::Vector3d b = point - vertices.col(i);


        // riduco a due equazioni (significative)
        // Eigen::Matrix2d Coeffs;
        // Eigen::Vector2d b_values;
        // if (abs(A.block(0, 0, 2, 2).determinant()) >= 5*numeric_limits<double>::epsilon()) {
        //     Coeffs = A.block(0, 0, 2, 2);
        //     b_values << b(0), b(1);
        // } else if (abs(A.block(1, 0, 2, 2).determinant()) >= 5*numeric_limits<double>::epsilon()) {
        //     Coeffs = A.block(1, 0, 2, 2);
        //     b_values << b(1), b(2);
        // } else {
        //     Coeffs << line(0), lato(0), line(2), lato(2); // controllare che funzioni
        //     b_values << b(0), b(2);
        // }
        //Eigen::Vector2d parameters = Coeffs.lu().solve(b_values);


        Eigen::Vector2d parameters = A.colPivHouseholderQr().solve(b);
        //se il parametro alfa è compreso tra 0 e 1, sarà verificata la combinazione convessa e troviamo il punto di intersezione
        if (0<=parameters(0) && parameters(0)<1) {
            Eigen::Vector3d point = vertices.col(i) + parameters(0)*lato;
            result.push_back(point);
        }
    }
    return result;
}


void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {
    /*stiamo lavorando su una frattura, passando in input un'altra frattura valutiamo se ci sia intersezione fra di esse e cataloghiamo il tipo di traccia ottenuta.
    Tale processo avviene su più step:
    -calcolo retta intersezione tra i piani contenenti le fratture
    -calcolo punti di intersezione tra la retta e i lati di entrambe le fratture
    -a seconda del numero totale di punti di intersezioni ottenuti, e alle rispettive distanze tra di essi, cataloghiamo il tipo di traccia

    Infine, tramite referenza, scriviamo i risultati sulla mesh delle tracce*/

    array<Eigen::Vector3d,2> trace_verteces;
    unsigned int trace_id = mesh.traces_id.size();

    //costruzione e risoluzione del sistema lineare per determinare la retta di intersezione tra i due piani contenenti le fratture
    //Per la risoluzione del sistema, poichè A è matrice quadrata; viene utilizzato il metodo PA=LU con costo computaziomale (n^3)/3
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

    //calcolo delle intersezioni tra i lati di entrambe le fratture e la retta ottenuta
    vector<Eigen::Vector3d> punti_1 = calculateIntersectionsPoints(t, p);
    vector<Eigen::Vector3d> punti_2 = other.calculateIntersectionsPoints(t, p);

    //definiamo il tipo di traccia a seconda del numero di punti intersezione ottenuti
    if (punti_1.size()+punti_2.size()!=4) {
        //se sono più di 4 o esattamente uno, è stato commesso un errore
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

    vector<Eigen::Vector3d> punti_distinti = Utils::calculateDistinctPoints(punti_1, punti_2); //vettore con tutti i punti di intersezione NON coincidenti

    if (punti_distinti.size() == 2) {
        //se sono esattamente 2, allora tali punti coincideranno (a due a due) su entrambe le fratture e la traccia sarà passante per entrambe
        array<Eigen::Vector3d, 2> punti_distinti_array; //punti di intersezione tra le due fratture
        copy_n(make_move_iterator(punti_distinti.begin()), 2, punti_distinti_array.begin());
        array<unsigned int, 2> fractures_id = {id, other.id};
        mesh.addTrace(trace_id, punti_distinti_array, fractures_id);
        passant_traces.push_back(id); //aggiungo l'id della traccia nel vettore delle tracce passanti (per la frattura su cui si sta lavorando)
        other.passant_traces.push_back(id); //aggiungo l'id della traccia nel vettore delle tracce passanti (per la frattura passata in input)
        return;
    }
    if (punti_distinti.size() == 3) {
        //se sono esattamente 3, allora due punti di intersezione in due fratture diverse, risulteranno coincidenti.
        //La traccia risulterà dunque passante per una frattura e interna per l'altra

        /*cerco quali siano i punti di intersezione coincidenti:
        per ogni punto di intersezione trovato, testo quale appartiene ad entrambe le fratture (aumentando il counter quando appartiene ad una frattura) */
        Eigen::Vector3d doppione; //coordinate del punto di intersezione che coincide
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

        //cerco quale sia la distanza minore tra il vertice in comune e l'altro punto di intersezione:
        //la traccia sarà passante per la frattura con distanza minore, e interna per l'altra
        double min_length = numeric_limits<double>::max();
        array<Eigen::Vector3d, 2> punti_distinti_array;
        for (Eigen::Vector3d& i: punti_distinti) {
            if (i == doppione) { //escludiamo la distanza vertice in comune/vertice in comune che vale 0 e sarà sempre la minore
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

        //verifico la distanza minore a quale frattura appartenga, e salvo di conseguenza la traccia come passante o interna
        //NB! sarà o passante nella frattura studiata e interna in quella passata in input, o viceversa
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



PolygonalMesh generatePolygonalMesh() {
    PolygonalMesh output;
    return output;
}
