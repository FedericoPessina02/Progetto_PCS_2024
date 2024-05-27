#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include "vector"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    // int partition_dimension = 2;
    // array<double, 6> domain_borders; //
    // vector<Fracture> fractures = Utils::fractureInput("./DFN/FR50_data.txt", domain_borders); //elenco fratture
    // map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, partition_dimension);
    // TracesMesh mesh;
    // Algorithms::cutTraces(id_to_fractures, mesh);
    // Utils::Stampa1("results1.csv",mesh);
    // Algorithms::ordinaFract(id_to_fractures, mesh,"results2.csv");
    // vector<PolygonalMesh> polygons = Algorithms::cutPolygonalMesh(id_to_fractures, mesh);

    ::testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}

// /*passant traces for each fracture;  internal traces;*/

// /*controllo correttezza traccia e correttezza tracce passanti e interne dati due poligoni*/
//perfeziona cosa dei poligoni messi in senso antiorario qualunque (basta che rispettino il senso di partenza)
TEST(zero_caso_interna_passante,norma_della_differenza_tra_punti_intersezione_corretti_e_calcolati) //caso - * * -
{
    unsigned int id1 = 26;
    unsigned int id2 = 9;
    unsigned int num_vertices1 = 4;
    unsigned int num_vertices2 = 3;
    TracesMesh mesh;
    MatrixXd vertices1;
    vertices1.resize(3,num_vertices1);
    //A=(-2,4,0)  B=(1,4,0)  C=(1,2,0)  D=(-2,2,0)
    vertices1 << -2,1,1,-2, 4,4,2,2, 0,0,0,0;
    MatrixXd vertices2;
    vertices2.resize(3,num_vertices2);
    //E=(0,4,1)  F=(0,2,1)  G=(0,3,-1)
    vertices2 << 0,0,0, 4,2,3, 1,1,-1;
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    Fracture h2 = Fracture(id2,num_vertices2,vertices2);
    h1.generateTrace(h2, mesh);
    Vector3d p1;
    Vector3d p2;
    Vector3d p1_real;
    p1_real << 0,2.5,0; /*sappiamo dai conti che esistiono e sono loro */
    Vector3d p2_real;
    p2_real << 0,3.5,0;
    for (unsigned int k = 0; k< 3; k++){
        p1[k] = mesh.traces_vertices[0][0][k];
        p2[k] = mesh.traces_vertices[0][1][k];
    }
    /* test per verificare correzza traccia */
    EXPECT_TRUE(((p1-p1_real).norm() + (p2-p2_real).norm() < 10 * numeric_limits<double>::epsilon())||((p1-p2_real).norm() + (p2-p1_real).norm() < 10 * numeric_limits<double>::epsilon()));
    /* test per verificare correttezza traccia passante o meno per le due date fratture */
    /* mi aspetto interna per h1 e passante per h2 */
    EXPECT_EQ(0,h1.passant_traces.size());
    EXPECT_EQ(1,h2.passant_traces.size());
    EXPECT_EQ(1,h1.internal_traces.size());
    EXPECT_EQ(0,h2.internal_traces.size());}

TEST(zero_traccia_inesistente,numero_di_tracce_salvate) //caso - - * *
{
    unsigned int id1 = 0;
    unsigned int id2 = 1;
    unsigned int num_vertices1 = 4;
    unsigned int num_vertices2 = 4;
    TracesMesh mesh;
    MatrixXd vertices1;
    vertices1.resize(3,num_vertices1);
    //A = 0 0 0;   B = 0,2,0;  C = 0,2,2;  D = 0 0 2;
    vertices1 << 0,0,0,0,0,2,2,0,0,0,2,2;
    MatrixXd vertices2;
    vertices2.resize(3,num_vertices2);
    //E=(4,0,1)  F=(4,2,1)  G=(3,1,1) H=(3,0,1)
    vertices2 << 4,4,3,3, 0,2,1,0, 1,1,1,1;
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    Fracture h2 = Fracture(id2,num_vertices2,vertices2);
    h1.generateTrace(h2, mesh);
    EXPECT_EQ(mesh.traces_id.size(),0); //non ho salvato tracce
    /*Osservazione: se una matrice con dimensione definita non viene riempita mi da come errore "more than 4 points of intersections */
}

TEST(zero_caso_interna_interna,norma_della_differenza_tra_punti_intersezione_corretti_e_calcolati) //caso - * - *
{
    unsigned int id1 = 110;
    unsigned int id2 = 11;
    unsigned int num_vertices1 = 4;
    unsigned int num_vertices2 = 4;
    TracesMesh mesh;
    MatrixXd vertices1;
    vertices1.resize(3,num_vertices1);
    //A=(-2,4,0)  B=(1,4,0)  C=(1,2,0)  D=(-2,2,0)
    vertices1 << -2,1,1,-2, 4,4,2,2, 0,0,0,0;
    MatrixXd vertices2;
    vertices2.resize(3,num_vertices2);
    //E=(0,3,-1)  F=(0,3,1)  G=(10,3,1)  H=(10,3,-1)
    vertices2 << 0,0,10,10, 3,3,3,3, -1,1,1,-1;
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    Fracture h2 = Fracture(id2,num_vertices2,vertices2);
    h1.generateTrace(h2, mesh);
    Vector3d p1;
    Vector3d p2;
    Vector3d p1_real;
    p1_real << 0,3,0; /* traccia calcolata by hand*/
    Vector3d p2_real;
    p2_real << 1,3,0;
    for (unsigned int k = 0; k< 3; k++){
        p1[k] = mesh.traces_vertices[0][0][k];
        p2[k] = mesh.traces_vertices[0][1][k];
    }
    /* test per verificare correttezza vertici traccia */
    EXPECT_TRUE(((p1-p1_real).norm() + (p2-p2_real).norm() < 10 * numeric_limits<double>::epsilon())||((p1-p2_real).norm() + (p2-p1_real).norm() < 10 * numeric_limits<double>::epsilon()));
    /* test per verificare correttezza traccia passante o meno per le due date fratture */
    /* mi aspetto sia interna sia per h1 che per h2 */
    EXPECT_EQ(0,h1.passant_traces.size());
    EXPECT_EQ(0,h2.passant_traces.size());
    EXPECT_EQ(1,h1.internal_traces.size());
    EXPECT_EQ(1,h2.internal_traces.size());
}

TEST(poligoni_tagliati_by_hand,poligoni_tagliati_da_cutPolygonalMesh)
{
    unsigned int id1 = 0;
    unsigned int num_vertices1 = 4;
    PolygonalMesh mesh;
    MatrixXd vertices1;
    vertices1.resize(3,num_vertices1);
    //punti poligono:
    Eigen::Vector3d A={-2,4,0};
    Eigen::Vector3d B={1,4,0};
    Eigen::Vector3d C={1,2,0};
    Eigen::Vector3d D={-2,2,0};
    vertices1 << -2,1,1,-2, 4,4,2,2, 0,0,0,0;
    //punti traccia:
    Eigen::Vector3d E={0,4,0};
    Eigen::Vector3d F={-1,2,0};
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    //ids:   A = 0;  B = 1;  C = 2;  D = 3;  E = 4;  F = 5;

    mesh.FractureId = id1;
    for (unsigned int i = 0; i < num_vertices1; i++) {
        mesh.CoordinateCell0Ds.push_back(vertices1.col(i));
        mesh.IdCell0Ds.push_back(i);
        if (i<num_vertices1-1) {
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
    unsigned int segment_point_1 = mesh.addPoint(E);
    unsigned int segment_point_2 = mesh.addPoint(F);
    array<unsigned int,2> segment = {segment_point_1,segment_point_2};
    array<unsigned int,2> intersection_starters = {0,2};
    Algorithms::cutPolygonBySegment(h1,mesh,id1,segment,intersection_starters); //ora vogliamo vedere se ha generato i due poligoni che ci aspettiamo
    //vogliamo scorrere sui poligoni e vedere i vertici di tali poligoni
    vector<unsigned int> vertici_giusti = {0,4,5,3,4,1,2,5};
    vector<unsigned int> vertici_da_algoritmo;
    vertici_da_algoritmo.reserve(8);
    unsigned int j = 0;
    for (unsigned int& polygon_id: mesh.activatedPolygons) {
        // cout << polygon_id << " ziopera" <<endl;
        // cout << endl;
        for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygon_id].size(); i ++){
            //cout << vertices1(i) << endl;
            //cout << mesh.VerticesCell2Ds[polygon_id][i] << endl;
            vertici_da_algoritmo.push_back(mesh.VerticesCell2Ds[polygon_id][i]);
            j += 1;
            // cout << j << endl;
        }
        // cout << endl;
    }
    if (vertici_da_algoritmo.size() != vertici_giusti.size()){
        ASSERT_EQ(0,1); //numero di vertici totali differente
    }else{
        ASSERT_TRUE(vertici_da_algoritmo == vertici_giusti);
    }
    //vertici del primo e secondo poligono rispettivamente coincidenti !
}

//Il test qua sotto non funziona, mi stampa dei poligoni che hanno piu volte lo stesso vertice -> DA MODIFICARE CUTPOLYGONBYSEGMENT
TEST(poligoni_tagliati_by_hand_2,poligoni_tagliati_da_cutPolygonalMesh) //la traccia ha come estremi i vertici del poligono
{
    //void cutPolygonBySegment(Fracture& fracture, PolygonalMesh& mesh, unsigned int polygonId, array<unsigned int,2> segment, array<unsigned int,2> intersection_starters) {
    unsigned int id1 = 0;
    unsigned int num_vertices1 = 5;
    PolygonalMesh mesh;
    MatrixXd vertices1;
    vertices1.resize(3,num_vertices1);
    //punti poligono:
    Eigen::Vector3d A={1,3,1};
    Eigen::Vector3d B={2,2,1};
    Eigen::Vector3d C={0,0,1};
    Eigen::Vector3d D={-1,1,1};
    Eigen::Vector3d E={-1,2,1};
    vertices1 << 1,2,0,-1,-1, 3,2,0,1,2, 1,1,1,1,1;
    //punti traccia:
    Eigen::Vector3d F={1,3,1};
    Eigen::Vector3d G={0,0,1};
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    //ids:   A = 0;  B = 1;  C = 2;  D = 3;  E = 4;  (F e G sono A e C)

    //in questo for costruisce piu punti (forse) controllare
    mesh.FractureId = id1;
    for (unsigned int i = 0; i < num_vertices1; i++) {
        mesh.CoordinateCell0Ds.push_back(vertices1.col(i));
        mesh.IdCell0Ds.push_back(i);
        if (i<num_vertices1-1) {
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
    // unsigned int segment_point_1 = mesh.addPoint(F);
    // unsigned int segment_point_2 = mesh.addPoint(G); sono 0 e 2 rispettivamente
    //cout << segment_point_1 << "   "  << segment_point_2 <<endl;
    array<unsigned int,2> segment = {0,2};
    array<unsigned int,2> intersection_starters = {0,2};
    Algorithms::cutPolygonBySegment(h1,mesh,id1,segment,intersection_starters); //ora vogliamo vedere se ha generato i due poligoni che ci aspettiamo
    //vogliamo scorrere sui poligoni e vedere i vertici di tali poligoni
    vector<unsigned int> vertici_giusti_1 = {0,1,2,2,3,4,0};
    vector<unsigned int> vertici_giusti_2 = {2,3,4,0,0,1,2}; //entrambi gli ordinamenti sono corretti
    vector<unsigned int> vertici_da_algoritmo;
    vertici_da_algoritmo.reserve(7);
    unsigned int j = 0;
    for (unsigned int& polygon_id: mesh.activatedPolygons) {
        for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygon_id].size(); i ++){
            vertici_da_algoritmo.push_back(mesh.VerticesCell2Ds[polygon_id][i]);
            j += 1;
            // cout << j << ":  ";
            // cout << mesh.VerticesCell2Ds[polygon_id][i] << endl;
        }
        // cout << endl;
    }
    if (vertici_da_algoritmo.size() != vertici_giusti_1.size()){
        ASSERT_EQ(0,1); //numero di vertici totali differente
    }else{
        ASSERT_TRUE(vertici_da_algoritmo == vertici_giusti_1 || vertici_da_algoritmo == vertici_giusti_2);
    }
   //vertici del primo e secondo poligono rispettivamente coincidenti !
}

//Test su cutMeshBySegment che ragiona solo con tracce passanti su una mesh (poligonale)
//fornisco un segmento che taglia due poligoni e mi aspetto come risultato 4 poligoni con i punti che li descrivono in senso antiorario
TEST(poligoni_corretti_e_ordinati_in_senso_antiorario,poligoni_forniti_da_algoritmo){
    //void Fracture::cutMeshBySegment(PolygonalMesh& mesh, Eigen::Vector3d a, Eigen::Vector3d b)
    unsigned int id0 = 0; //della frattura
    unsigned int id1 = 1; //dei sottopoligoni
    unsigned int id2 = 2;
    unsigned int num_vertices0 = 4;
    unsigned int num_vertices1 = 4;
    unsigned int num_vertices2 = 4;
    PolygonalMesh mesh;
    MatrixXd vertices0;
    MatrixXd vertices1;
    MatrixXd vertices2;
    Eigen::Vector3d A={0,0,0}; //0
    Eigen::Vector3d B={3,0,0}; //4
    Eigen::Vector3d C={6,0,0}; //1
    Eigen::Vector3d D={6,3,0}; //2
    Eigen::Vector3d E={3,3,0}; //5
    Eigen::Vector3d F={0,3,0}; //3
    // A C D F lati poligonone
    vertices0.resize(3,num_vertices0);
    vertices0 << 0,6,6,0, 0,0,3,3, 0,0,0,0;

    Fracture h0 = Fracture(id0,num_vertices0,vertices0);
    //la frattura è quella data da h0
    mesh.FractureId = id0;
    for (unsigned int i = 0; i < num_vertices0; i++) {
        mesh.CoordinateCell0Ds.push_back(vertices0.col(i));
        mesh.IdCell0Ds.push_back(i);
        if (i<num_vertices0-1) {
            mesh.IdCell1Ds.push_back(mesh.IdCell1Ds.size());
            mesh.VerticesCell1Ds.push_back(array<unsigned int,2> {i, i+1});
        } else {
            mesh.IdCell1Ds.push_back(mesh.IdCell1Ds.size());
            mesh.VerticesCell1Ds.push_back(array<unsigned int,2> {i, 0});
        }
    }
    mesh.IdCell2Ds.push_back(0);
    mesh.VerticesCell2Ds.push_back(mesh.IdCell0Ds);
    mesh.EdgesCell2Ds.push_back(mesh.IdCell1Ds);

    //inserimento punti necessari nella mesh
    mesh.addPoint(B);
    mesh.addPoint(E);

    //inserimento lati per poligoni 1 e 2 nella mesh
    // 1
    mesh.addEdge(0,4); mesh.addEdge(4,5); mesh.addEdge(5,3); //lati di coordinate 4 5 6
    //2
    mesh.addEdge(4,1); mesh.addEdge(2,5); mesh.addEdge(5,4); //lati di coordinate 7 8 9

    mesh.IdCell2Ds.push_back(id1);
    mesh.IdCell2Ds.push_back(id2);
    //aggiungo punti a IdCell2Ds e anche lati
    vector<unsigned int> pol1_vertices = {0,4,5,3}; // A B E F
    vector<unsigned int> pol2_vertices = {4,1,2,5}; // B C D E
    vector<unsigned int> pol1_edges = {4,5,6,3}; // AB BE EF FA
    vector<unsigned int> pol2_edges = {7,1,8,9}; // BC CD DE EB
    mesh.VerticesCell2Ds.push_back(pol1_vertices);
    mesh.VerticesCell2Ds.push_back(pol2_vertices);
    mesh.EdgesCell2Ds.push_back(pol1_edges);
    mesh.EdgesCell2Ds.push_back(pol2_edges);
    mesh.activatedPolygons.push_back(id1);
    mesh.activatedPolygons.push_back(id2);

    Eigen::Vector3d a={6,2,0};
    Eigen::Vector3d b={0,1,0};

    h0.cutMeshBySegment(mesh,a,b);
    cout << mesh.activatedPolygons.size() << endl;
    for (unsigned int i = 0; i < mesh.activatedPolygons.size(); i++){
        cout << mesh.activatedPolygons[i] << ": ";
        for (unsigned int j = 0; j < mesh.VerticesCell2Ds[i].size(); j++){
            cout << mesh.VerticesCell2Ds[i][j] << "   ";
        }
        cout << endl;
    }
}

//Manca test su generatePolygonalMesh che fa tutto quindi prendi un esempio di poligono con 2 tracce e dentro e controlla tutti i tagli corretti

/* */

