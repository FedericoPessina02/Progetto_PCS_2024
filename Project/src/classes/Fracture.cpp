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

void Fracture::generateTrace(Fracture& other, TracesMesh& mesh) {
    /* stiamo lavorando su una frattura, passando in input un'altra frattura valutiamo se ci sia intersezione fra di esse e cataloghiamo il tipo di traccia ottenuta.
    Tale processo avviene su più step:
    -calcolo retta intersezione tra i piani contenenti le fratture
    -calcolo punti di intersezione tra la retta e i lati di entrambe le fratture
    -a seconda del numero totale di punti di intersezioni ottenuti, e alle rispettive distanze tra di essi, cataloghiamo il tipo di traccia

    Infine, tramite referenza, scriviamo i risultati sulla mesh delle tracce */

    array<Eigen::Vector3d,2> trace_verteces;
    unsigned int trace_id = mesh.traces_id.size();

    /* costruzione e risoluzione del sistema lineare per determinare la retta di intersezione tra i due piani contenenti le fratture
    Per la risoluzione del sistema, poichè A è matrice quadrata; viene utilizzato il metodo PA=LU con costo computaziomale (n^3)/3 */
    Eigen::Vector3d t = normal.cross(other.normal);
    Eigen::Matrix3d A;
    A.row(0) = normal;
    A.row(1) = other.normal;
    A.row(2) = t;
    if (abs(A.determinant()) < 5*numeric_limits<double>::epsilon()) {
        // la matrice associata al sistema è singolare, quindi il sistema non è determinato (ossia le due fratture sono parallele nello spazio)
        return;
    }
    Eigen::Vector3d b(plane_d, other.plane_d, 0); // <-- vettore colonna contenente i termini noti delle equazioni planari associate
    Eigen::Vector3d p = A.lu().solve(b);

    // una volta ricavata l'equazione parametrica della retta lungo cui corre la traccia calcolo i punti di intersezione
    // ogni frattura restituisce una lista di punti di intersezione trovati su di lei
    vector<Eigen::Vector3d> punti_1 = Algorithms::calculateIntersectionsPoints(*this, t, p);
    vector<Eigen::Vector3d> punti_2 = Algorithms::calculateIntersectionsPoints(other, t, p);

    // definiamo il tipo di traccia a seconda del numero di punti intersezione ottenuti
    // escludo il caso in cui vi sia un singolo punto di intersezione (traccia tangente a un poligono)
    // se vi sono più di 4 punti segnalo l'errore (i poligoni sono convessi non può verificarsi questa condizione)
    if (punti_1.size()+punti_2.size()!=4) {
        if (punti_1.size()+punti_2.size()>4) {
            cerr << "More than 4 points of intersections";
            return;
        }
        return;
    }

    // genero un vettore contenente l'unione distinta dei punti di intersezione
    // la sua lunghezza mi permette di discriminare quale algoritmo utilizzare
    // il caso generale con 4 punti di intersezione andrebbe bene sempre, ma questa verifica mi permette di ridurre se possibile i calcoli
    //     sia in termini di controlli fatti che di calcoli svolti
    vector<Eigen::Vector3d> punti_distinti = Utils::calculateDistinctPoints(punti_1, punti_2); //vettore con tutti i punti di intersezione NON coincidenti

    if (punti_distinti.size() == 2) {
        // se sono esattamente 2, allora tali punti coincideranno (a due a due) su entrambe le fratture e la traccia sarà passante per entrambe
        // in questo caso ho la garanzia che la traccia esiste e non vi possono essere casi degeneri
        array<Eigen::Vector3d, 2> punti_distinti_array; //punti di intersezione tra le due fratture
        copy_n(make_move_iterator(punti_distinti.begin()), 2, punti_distinti_array.begin());
        array<unsigned int, 2> fractures_id = {id, other.id};
        mesh.addTrace(trace_id, punti_distinti_array, fractures_id);
        passant_traces.push_back(trace_id); //aggiungo l'id della traccia nel vettore delle tracce passanti (per la frattura su cui si sta lavorando)
        other.passant_traces.push_back(trace_id); //aggiungo l'id della traccia nel vettore delle tracce passanti (per la frattura passata in input)
        return;
    }
    if (punti_distinti.size() == 3) {
        // se sono esattamente 3, allora due punti di intersezione in due fratture diverse, risulteranno coincidenti.
        // la traccia risulterà dunque passante per una frattura e interna per l'altra
        // come nel caso prima ho la garanzia che la traccia esiste e non vi possono essere casi degeneri
        //      (i casi di tracce tangenti ai poligoni li ho esclusi prima)

        // cerco quali siano i punti di intersezione coincidenti:
        // per ogni punto di intersezione trovato, testo quale appartiene ad entrambe le fratture
        //      (aumentando il counter quando appartiene ad una frattura)
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

        // cerco quale sia la distanza minore tra il vertice in comune e l'altro punto di intersezione:
        // la traccia sarà passante per la frattura con distanza minore, e interna per l'altra
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
        //NB sarà o passante nella frattura studiata e interna in quella passata in input, o viceversa
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

    // il caso generale a 4 punti richiede di ordinare i punti e di studiare l'appartenenza dei vari punti ai poligoni
    // come prima cosa individuo un punto estremale (salvato in estremo) e memorizzo la referenza a quale poligoni appartiene
    //     so per certo che appartiene a un solo poligono perché altrimenti non starei usando questo algoritmo
    vector<Eigen::Vector3d>* punti_ptr = &punti_1;
    double distanza = 0;
    Eigen::Vector3d estremo;
    for (unsigned int i = 0; i < punti_1.size(); i++) {
        for (unsigned int j = 0; j < punti_2.size(); j++) {
            if ((punti_1[i]-punti_2[j]).squaredNorm() > distanza) {
                estremo = punti_1[i];
                distanza = (punti_1[i]-punti_2[j]).squaredNorm();
            }
        }
    }
    if ((punti_2[0]-punti_2[1]).squaredNorm() > distanza) {
        estremo = punti_2[0];
        punti_ptr = &punti_2;
    }

    // a partire il punto estremale salvo in punto_vicini il punto a lui più vicino
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

    // controllo il caso degenere di due fratture separate con segmento di intersezione esterno ad entrambi i poligoni
    // se punto_vicini appartiene alla coppia di punti a cui appartiene anche il punto estremale che sto considerando ricado nel caso degenere:
    //      le due fratture non si stanno realmente intersecando
    if ((punti_ptr->operator[](0)-punto_vicini).squaredNorm() < 8*numeric_limits<double>::epsilon()) {
        return;
    }
    if ((punti_ptr->operator[](1)-punto_vicini).squaredNorm() < 8*numeric_limits<double>::epsilon()) {
        return;
    }

    // cerco il secondo punto più vicino al punto estremale
    // NB esso non corrisponde necessariamente al punto più vicino a punto_vicino, è solo il secondo punto più vicino al punto estremale scelto
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

    // la coppia {punto_vicini, punto_meno_vicini} rappresenta gli estremi della traccia trovata
    // procedo a salvare la traccia
    vector<Eigen::Vector3d> punti_traccia_vector = {punto_vicini, punto_meno_vicini};
    array<Eigen::Vector3d, 2> punti_traccia_array = {punto_vicini, punto_meno_vicini};
    array<unsigned int, 2> fractures_id = {id, other.id};
    mesh.addTrace(trace_id, punti_traccia_array, fractures_id);

    // essendo questo il caso generale non posso concludere a priori la natura della traccia
    // quindi per ogni frattura contronto gli estremi della traccia con i punti di intersezione da lei generati
    // se coincidono la traccia è passante, altrimenti la traccia è interna
    if (Utils::compareSegments(punti_traccia_vector, punti_1)) {
        passant_traces.push_back(trace_id);
    } else {
        internal_traces.push_back(trace_id);
    }
    if (Utils::compareSegments(punti_traccia_vector, punti_2)) {
        other.passant_traces.push_back(trace_id);
    } else {
        other.internal_traces.push_back(trace_id);
    }
}

PolygonalMesh Fracture::generatePolygonalMesh(TracesMesh& traces) {
    /* il metodo gestisce la generazione della mesh poligonale a partire dai dati disponibili:
    la struttura della frattura (che costituirà il primo poligono della mesh) e l'elenco delle fratture
    esse sono divise in passanti e interne (si assumono già ordinate per lunghezza decrescente) */

    // come prima cosa viene creata la mesh poligonale "vuota"
    // essa contiene un solo poligono, il quale ricalca la struttura della frattura
    PolygonalMesh mesh;
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
    // per ipotesi gli estremi delle tracce sono punti appartenenti a dei segmenti della mesh poligonale, di conseguenza
    //      posso usarli direttamente come estremi di taglio
    for (unsigned int& trace_id: passant_traces) {
        cutMeshBySegment(mesh, traces.traces_vertices[trace_id][0], traces.traces_vertices[trace_id][1]);
    }

    // taglio per tracce non passanti
    // è necessario prolungare la traccia fino a giungere a un segmento che rispetti le ipotesi necessarie per il taglio
    // estendo la traccia lungo la sua retta e prendo tutti i punti di intersezione, procedo poi a cercare i punti esterni alla traccia
    //     aventi rispettivamente il coefficiente <= 0 più grande e il coefficiente >= 1 più piccolo
    //     essi corrispondono agli estremi della traccia passante su cui lanciare il taglio
    for (unsigned int& trace_id: internal_traces) {
        // creo la struttura dati che conterrà le soluzioni dei sistemi lineari già calcolati
        // unordered_map è tabella hash e come tale ha un tempo di accesso diretto O(1)/O(N) (nel caso peggiore di sole collisioni)
        // la map tradizionale è implementata con una struttura simile ad un albero AVL, quindi ha tempo di accesso O(logN)
        unordered_map<unsigned int, array<double,2>> cached_coeffs;
        // trovo l'equazione parametrica della retta che estende la traccia
        Eigen::Vector3d direction = traces.traces_vertices[trace_id][1] - traces.traces_vertices[trace_id][0];
        Eigen::Vector3d application_point = traces.traces_vertices[trace_id][0];

        vector<Eigen::Vector3d> standard_intersection_points; // punti di intersezione lungo l'estensione della traccia
        vector<double> combination_coeffs; // coefficienti dell'equazione parametrica della retta associati ai rispetti punti
        for (unsigned int& polygonId: mesh.activatedPolygons) {
            for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygonId].size(); i++) {
                // scorro i lati del poligono per cercare i punti di intersezione
                unsigned int vertex_a_id = mesh.VerticesCell2Ds[polygonId][i];
                unsigned int vertex_b_id;
                Eigen::Vector3d a = mesh.CoordinateCell0Ds[vertex_a_id];
                Eigen::Vector3d b;
                if (i < mesh.VerticesCell2Ds[polygonId].size() - 1) {
                    vertex_b_id = mesh.VerticesCell2Ds[polygonId][i+1];
                    b = mesh.CoordinateCell0Ds[vertex_b_id];
                } else {
                    vertex_b_id = mesh.VerticesCell2Ds[polygonId][0];
                    b = mesh.CoordinateCell0Ds[vertex_b_id];
                }
                // trovo l'equazione parametrica associata al lato del poligono
                Eigen::Vector3d edge_direction = b-a;
                Eigen::Vector3d edge_application = a;
                if (edge_direction.cross(direction).norm() < 5*numeric_limits<double>::epsilon()) {
                    // il lato è parallelo al segmento, non ha senso cercare una possibile intersezione
                    continue;
                }
                // imposto il sistema (sovradeterminato) per cercare l'intersezione
                Eigen::MatrixXd A;
                A.resize(3,2);
                A.col(0) = edge_direction;
                A.col(1) = -1*direction;
                Eigen::Vector3d coef = application_point - edge_application;
                Eigen::Vector2d parameters = A.colPivHouseholderQr().solve(coef);
                Eigen::Vector3d point = a + parameters(0)*edge_direction;
                if (-5*numeric_limits<double>::epsilon()<=parameters(0) && parameters(0)<1+5*numeric_limits<double>::epsilon()) {
                    // se la condizione è rispettata vuol dire che l'intersezione cade dentro il lato e non nella sua estesa
                    // quindi posso considerarlo come intersezione (insieme all'altro coefficiente, associato all'equazione parametrica
                    // della retta su cui giace la traccia)
                    unsigned int edge_id = mesh.findEdge(vertex_a_id, vertex_b_id);
                    cached_coeffs[edge_id] = array{parameters(0), parameters(1)}; // memorizzo la soluzione del sistema lineare
                    standard_intersection_points.push_back(point);
                    combination_coeffs.push_back(parameters(1));
                }
            }
        }

        // distingue tre casi --> coefficienti tutti dello stesso segno o misti
        // questa distinzione è dovuta al fatto che non posso a priori distinguere l'orientazione della direttrice
        // probabilmente si può fare diversamente, ma così funziona
        array<Eigen::Vector3d, 2> extended_trace_borders;
        unsigned int strict_positives = 0;
        unsigned int strict_negatives = 0;
        for (double& coeff: combination_coeffs) {
            if (coeff > 0) {
                strict_positives ++;
            } else if (coeff < 0) {
                strict_negatives ++;
            }
        }

        if (strict_positives != 0 && strict_negatives != 0) {
            // se alcuni coefficienti sono negativi e altri sono positivi mi ritrovo nel caso standard
            // trovo i due punti estremali che corrispondono alla minore estensione passante possibile della traccia
            // per farlo cerco il punto con coefficiente negativo più alto (attenzione che può essere 0 se quell'estremo è effettivamente su un lato)
            // e quello con il coefficiente positivo e maggiore o uguale a uno più alto (stesso discorso di prima ma con le disuguaglianze capovolte)
            double negative_coeff = numeric_limits<double>::lowest();
            double positive_coeff = numeric_limits<double>::max();
            for (unsigned int i = 0; i<standard_intersection_points.size(); i++) {
                if (combination_coeffs[i] <= 0) {
                    if (combination_coeffs[i] > negative_coeff) {
                        extended_trace_borders[0] = standard_intersection_points[i];
                        negative_coeff = combination_coeffs[i];
                    }
                } else if (combination_coeffs[i] >= 1) {
                    if (combination_coeffs[i] < positive_coeff) {
                        extended_trace_borders[1] = standard_intersection_points[i];
                        positive_coeff = combination_coeffs[i];
                    }
                }
            }
        } else {
            // se i segni sono tutti concordi vuol dire che la direttrice della retta parametrica punta al punto di applicazione
            // il punto associato al coefficiente 0 sarà sempre un estremo dell'estesa (si verifica manualmente)
            // per il secondo estremo cerco il punto che minimizza effettivamente la distanza senza però ricadere in intersezioni con la traccia interna
            double min_coeff = numeric_limits<double>::max();
            for (unsigned int i = 0; i<standard_intersection_points.size(); i++) {
                if (combination_coeffs[i] == 0) {
                    extended_trace_borders[0] = standard_intersection_points[i];
                } else if (abs(combination_coeffs[i]) < min_coeff && abs(combination_coeffs[i] >= 1)) {
                    extended_trace_borders[1] = standard_intersection_points[i];
                    min_coeff = abs(combination_coeffs[i]);
                }
            }
        }

        // lancio il taglio della mesh usando come estremi l'estesa della traccia (che a questo punto è passante)
        cutMeshBySegment(mesh, extended_trace_borders[0], extended_trace_borders[1]);
    }

    return mesh;
}

void Fracture::cutMeshBySegment(PolygonalMesh& mesh, Eigen::Vector3d& a, Eigen::Vector3d& b, const unordered_map<unsigned int, array<double,2>>& cached_coeffs) {
    // taglia la mesh intera lungo un segmento che ha come estremi due punti appartenenti a due lati dentro la mesh
    // cerco tutte le intersezioni e lancio il taglio su tutti i poligoni interessati usando come estremi di taglio il segmento della traccia che
    //      passa lungo quel poligono
    // come prima cosa trovo la retta su cui giace il segmento di taglio
    bool use_cached_coeffs = true; // se true il programma prova prima a vedere se le soluzioni dei sistemi lineari sono già state calcolate in precedenza
    Eigen::Vector3d direction = b-a;
    Eigen::Vector3d application_point = a;
    vector<unsigned int> to_be_modified_polygons = mesh.activatedPolygons; //id dei poligoni interessati dal taglio
    for (unsigned int& polygonId: to_be_modified_polygons) {
        // voglio memorizzare per ogni poligono le intersezioni dei suoi lati con il segmento e memorizzarmi anche tra quali vertici cadono
        //      per questa seconda informazione mi basta memorizzare i vertici che li precedono in senso antiorario (salvati in intersection_starters)
        vector<Eigen::Vector3d> intersection_points;
        vector<unsigned int> intersection_starters;
        // scorro i lati del poligono per cercare i punti di intersezione
        // il codice è analogo a quello contenuto in cutPolygonalMesh
        for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygonId].size(); i++) {
            unsigned int vertex_a_id = mesh.VerticesCell2Ds[polygonId][i];
            unsigned int vertex_b_id;
            Eigen::Vector3d a = mesh.CoordinateCell0Ds[vertex_a_id];
            Eigen::Vector3d b;
            if (i < mesh.VerticesCell2Ds[polygonId].size() - 1) {
                vertex_b_id = mesh.VerticesCell2Ds[polygonId][i+1];
                b = mesh.CoordinateCell0Ds[vertex_b_id];
            } else {
                vertex_b_id = mesh.VerticesCell2Ds[polygonId][0];
                b = mesh.CoordinateCell0Ds[vertex_b_id];
            }
            Eigen::Vector3d edge_direction = b-a;
            Eigen::Vector3d edge_application = a;
            if (edge_direction.cross(direction).norm() < 5*numeric_limits<double>::epsilon()) {
                continue; // il lato è parallelo al segmento non ha senso cercare un'intersezione
            }
            Eigen::Vector2d parameters;
            Eigen::MatrixXd A;

            unsigned int edge_id = mesh.findEdge(vertex_a_id, vertex_b_id);
            if (cached_coeffs.find(edge_id) != cached_coeffs.end() && use_cached_coeffs) {
                // se esiste la chiave associata al lato vuol dire che è già stato risolto il sistema lineare associato (bisogna essere sicuri che sia quello giusto ovviamente...)
                 // .at() è necessario perché cached_coeffs è dichiarato come const
                // cached_coeffs è dichiarato come const per permettere di assegnare un valore di default senza una reference (vedere Fracture.hpp)
                parameters(0) = cached_coeffs.at(edge_id)[0], cached_coeffs.at(edge_id)[1];
            } else {
                A.resize(3,2);
                A.col(0) = edge_direction;
                A.col(1) = -1*direction;
                Eigen::Vector3d coef = application_point - edge_application;
                parameters = A.colPivHouseholderQr().solve(coef);
            }
            if (-5*numeric_limits<double>::epsilon()<=parameters(0) && parameters(0)<1-5*numeric_limits<double>::epsilon()) {
                if (-5*numeric_limits<double>::epsilon()<=parameters(1) && parameters(1)<1+5*numeric_limits<double>::epsilon()) {
                    // oltre a verificare che l'intersezione sia propria (ossia che cada dentro il lato e non nella sua estesa)
                    // devo anche assicurarmi che si verifichi la stessa cosa per il segmento di taglio (ossia che vada a tagliare fuori)
                    // questo controllo è necessario perché se non ci fosse vanificherebbe tutto il lavoro fatto per minimizzare l'estensione delle tracce interne
                    Eigen::Vector3d point = a + parameters(0)*edge_direction;
                    intersection_points.push_back(point); // memorizzo il punto d'intersezione
                    intersection_starters.push_back(mesh.VerticesCell2Ds[polygonId][i]); // memorizzo il punto preesistente che lo precede in senso antiorario
                }
            }
        }
        if (intersection_points.size() != 2) {
            // non ho trovato abbastanza intersezioni, quindi il poligono non è interessato al taglio e posso skippare il taglio effettivo
            continue;
        }
        array<unsigned int, 2> segment; // conterrà gli id degli estremi del segmento di taglio
        // NB segment ora contiene un sottointervallo del segmento di taglio, ossia quel sottointervallo che effettivamente interessa il poligono in questione
        // memorizzo nella mesh gli estremi di taglio e mi salvo i loro id (se esistono già mi restituirà gli id esistenti)
        segment[0] = mesh.addPoint(intersection_points[0]);
        segment[1] = mesh.addPoint(intersection_points[1]);
        // creo ora un vettore contenente gli id di tutti i punti (ossia vado a inserire i punti nuovi originati dal taglio nelle posizioni giuste)
        // se il punto C originato dal taglio passa sul lato AB il risultato sarà A-C-B (la conoscenza degli intersection_starters mi preserva l'ordinamento in senso antiorario
        vector<unsigned int> total_points;
        for (unsigned int& vertex_id: mesh.VerticesCell2Ds[polygonId]) {
            total_points.push_back(vertex_id);
            if (vertex_id == intersection_starters[0] && vertex_id != segment[0] ) {
                total_points.push_back(segment[0]);
            } else if (vertex_id == intersection_starters[1] && vertex_id != segment[1] ) {
                total_points.push_back(segment[1]);
            }
        }
        // ora che ho ridotto il problema nella sua forma più semplice posso lanciare il taglio effettivo del poligono
        Algorithms::cutPolygonBySegment(*this, mesh, polygonId, total_points, segment);
    }

}
