#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;
using namespace Eigen;

TEST(calculateSphereTest,radius)
{   //NB. Il raggio calcolato è al quadrato
    unsigned int id = 0;
    unsigned int num_vertices = 4;
    Eigen::MatrixXd poli1;
    poli1.resize(3,4);
    poli1 << 0,1,1,0,0,0,1,1,0,0,0,0;
    Fracture frattura(id,num_vertices,poli1);
    frattura.calculateSphere();

    unsigned int id2 = 1;
    unsigned int num_vertices2 = 3;
    Eigen::MatrixXd poli2;
    poli2.resize(3,3);
    poli2 << 0,3,1.5,0,0,2.5980762113533,0,0,0;
    Fracture frattura2(id2,num_vertices2,poli2);
    frattura2.calculateSphere();


    EXPECT_TRUE(abs(frattura.radius-0.5)<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(frattura2.radius-3)<Utils::tol_coeff*numeric_limits<double>::epsilon());

}


TEST(calculateNormalVectorTest, normal)
{
    unsigned int id2 = 0;
    unsigned int num_vertices2 = 5;
    Eigen::MatrixXd poli2;
    poli2.resize(3,5);
    poli2 << 0,2,0,4,7,0,0,1,9,8,0,0,0,0,0;
    Fracture frattura2(id2,num_vertices2,poli2);
    frattura2.calculateNormalVector();

    EXPECT_TRUE(abs(0-frattura2.normal[0])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(0-frattura2.normal[1])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(2-frattura2.normal[2])<Utils::tol_coeff*numeric_limits<double>::epsilon());

    EXPECT_TRUE(abs(0-frattura2.plane_d)<Utils::tol_coeff*numeric_limits<double>::epsilon());

}


TEST(calculateIntersectionsPointsTest, intersections)
{
    unsigned int id = 0;
    unsigned int num_vertices = 4;
    Eigen::MatrixXd poli;
    poli.resize(3,4);
    poli << 0,2,8,1,0,-2,4,3,0,0,0,0;
    Fracture frattura(id,num_vertices,poli);
    Eigen::Vector3d line;
    line << 1,2,0;
    Eigen::Vector3d point;
    point << 0,-1,0;
    vector<Eigen::Vector3d> intersezioni= Algorithms::calculateIntersectionsPoints(frattura, line, point);

    unsigned int id1 = 1;
    unsigned int num_vertices1 = 4;
    Eigen::MatrixXd poli1;
    poli1.resize(3,4);
    poli1 << 0,1,1,0,0,0,1,1,0,0,0,0;
    Fracture frattura1(id1,num_vertices1,poli1);
    Eigen::Vector3d line1;
    line1 << 1,0,0;
    Eigen::Vector3d point1;
    point1 << 0,-1,0;
    vector<Eigen::Vector3d> intersezioni1= Algorithms::calculateIntersectionsPoints(frattura1, line1, point1);

    //caso con intersezioni
    EXPECT_TRUE(abs(0.333333333333333333333333333333333-intersezioni[0][0])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(-0.33333333333333333333333333333333-intersezioni[0][1])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(0-intersezioni[0][2])<Utils::tol_coeff*numeric_limits<double>::epsilon());

    EXPECT_TRUE(abs(2.0769230769231-intersezioni[1][0])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(3.1538461538462-intersezioni[1][1])<Utils::tol_coeff*numeric_limits<double>::epsilon());
    EXPECT_TRUE(abs(0-intersezioni[1][2])<Utils::tol_coeff*numeric_limits<double>::epsilon());

    //caso senza intersezioni
    ASSERT_EQ(0,intersezioni1.size());
}


TEST(assignPartitionTest, partitions)
{
    vector<Fracture> fractures;
    fractures.reserve(3);

    array<double, 6> domain_borders = {0, 10, 0, 10, 0, 10};

    const int partitions_number=2;

    //frattura nella prima partizione
    unsigned int id1 = 0;
    unsigned int num_vertices1 = 4;
    Eigen::MatrixXd poli1;
    poli1.resize(3,4);
    poli1 << 1,2,1,4,1,0,0,2,1,4,3,2;
    Fracture frattura1(id1,num_vertices1,poli1);
    fractures.push_back(frattura1);

    //frattura nell'ultima partizione
    unsigned int id2 = 1;
    unsigned int num_vertices2 = 4;
    Eigen::MatrixXd poli2;
    poli2.resize(3,4);
    poli2 << 7,8,9,9,8,8,6,9,9,8,7,7;
    Fracture frattura2(id2,num_vertices2,poli2);
    fractures.push_back(frattura2);

    //frattura a cavallo tra le partizioni
    unsigned int id3 = 2;
    unsigned int num_vertices3 = 4;
    Eigen::MatrixXd poli3;
    poli3.resize(3,4);
    poli3 << 4,6,5,3,4,6,6,2,4,6,7,1;
    Fracture frattura3(id3,num_vertices3,poli3);
    fractures.push_back(frattura3);

    map<int, vector<Fracture>> id_to_fractures=Algorithms::assignPartition(fractures, domain_borders, partitions_number);

    vector<Fracture> vettore_fratture_partizione;
    vettore_fratture_partizione=id_to_fractures[0];
    ASSERT_EQ(vettore_fratture_partizione.size(),1);
    ASSERT_EQ(vettore_fratture_partizione[0].id,frattura3.id);

    vettore_fratture_partizione=id_to_fractures[1];
    ASSERT_EQ(vettore_fratture_partizione.size(),1);
    ASSERT_EQ(vettore_fratture_partizione[0].id,frattura1.id);

    vettore_fratture_partizione=id_to_fractures[8];
    ASSERT_EQ(vettore_fratture_partizione.size(),1);
    ASSERT_EQ(vettore_fratture_partizione[0].id,frattura2.id);

}


TEST(compareSegmentsTest,comparison)
{
    vector<Eigen::Vector3d> vec1(2), vec2(2), vec3(2), vec4(2);
    vec1[0] << 5,8,4;
    vec1[1] << 1,1,1;

    vec2[0] << 5,8,4;
    vec2[1] << 1,1,1;

    vec3[0] << 1,1,1;
    vec3[1] << 5,8,4;

    vec4[0] << 2,3,1;
    vec4[1] << 5,8,4;

    ASSERT_EQ(true,Utils::compareSegments(vec1,vec2)); //entrambi i vertici coincidono
    ASSERT_EQ(true,Utils::compareSegments(vec1,vec3)); //vertici invertiti ma con stesse coordinate
    ASSERT_EQ(false,Utils::compareSegments(vec1,vec4)); //vertici distinti

}


TEST(calculateDistinctPointsTest,points)
{
    vector<Eigen::Vector3d> a(2), b(2), c(2), d(2);
    a[0] << 5,8,4;
    a[1] << 1,1,1;

    b[0] << 5,8,4;
    b[1] << 1,1,1;

    c[0] << 1,2,1;
    c[1] << 5,8,4;

    d[0] << 2,3,1;
    d[1] << 5,4,4;

    ASSERT_EQ(2,Utils::calculateDistinctPoints(a,b).size()); //2 vertici distinti
    ASSERT_EQ(3,Utils::calculateDistinctPoints(a,c).size()); //3 vertici distinti
    ASSERT_EQ(4,Utils::calculateDistinctPoints(a,d).size()); //4 vertici distinti
}

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
    EXPECT_TRUE(((p1-p1_real).norm() + (p2-p2_real).norm() < Utils::tol_coeff * numeric_limits<double>::epsilon())||((p1-p2_real).norm() + (p2-p1_real).norm() < Utils::tol_coeff * numeric_limits<double>::epsilon()));
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
    EXPECT_TRUE(((p1-p1_real).norm() + (p2-p2_real).norm() < Utils::tol_coeff * numeric_limits<double>::epsilon())||((p1-p2_real).norm() + (p2-p1_real).norm() < Utils::tol_coeff * numeric_limits<double>::epsilon()));
    /* test per verificare correttezza traccia passante o meno per le due date fratture */
    /* mi aspetto sia interna sia per h1 che per h2 */
    EXPECT_EQ(0,h1.passant_traces.size());
    EXPECT_EQ(0,h2.passant_traces.size());
    EXPECT_EQ(1,h1.internal_traces.size());
    EXPECT_EQ(1,h2.internal_traces.size());
}

TEST(poligoni_tagliati_by_hand,poligoni_tagliati_da_cutPolygonBySegment) //cutPolygonbysegment   FUNZIONA
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
    // Eigen::Vector3d E={1,4,0}; //CASO TRACCIA PASSANTE PER ANGOLI
    // Eigen::Vector3d F={-2,2,0};
    Fracture h1 = Fracture(id1,num_vertices1,vertices1);
    //ids:   A = 0;  B = 1;  C = 2;  D = 3;  E = 4;  F = 5;

    //aggiungo punti e vettori nella mesh
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
    vector<unsigned int> total_points = {0,4,1,2,5,3};
    //vector<unsigned int> total_points = {0,1,2,3}; //CASO TRACCIA PASSANTE PER ANGOLI
    Algorithms::cutPolygonBySegment(h1,mesh,id1,total_points,segment); //ora vogliamo vedere se ha generato i due poligoni che ci aspettiamo
    //vogliamo scorrere sui poligoni e vedere i vertici di tali poligoni
    vector<unsigned int> vertici_giusti = {0,4,5,3,4,1,2,5};
    //vector<unsigned int> vertici_giusti = {0,1,3,1,2,3}; //CASO TRACCIA PASSANTE PER ANGOLI
    vector<unsigned int> vertici_da_algoritmo;
    vertici_da_algoritmo.reserve(8);
    unsigned int j = 0;
    for (unsigned int& polygon_id: mesh.activatedPolygons) {
        cout << polygon_id << endl;
        for (unsigned int i = 0; i < mesh.VerticesCell2Ds[polygon_id].size(); i ++){
            vertici_da_algoritmo.push_back(mesh.VerticesCell2Ds[polygon_id][i]);
            cout << "  " << mesh.VerticesCell2Ds[polygon_id][i];
        }
         cout << endl;
    }
    if (vertici_da_algoritmo.size() != vertici_giusti.size()){
        ASSERT_EQ(0,1); //numero di vertici totali differente
    }else{
        ASSERT_TRUE(vertici_da_algoritmo == vertici_giusti);
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
    Eigen::Vector3d b = {0,1,0};
    //Eigen::Vector3d a = {3,3,0}; //CASO TRACCIA PASSANTE PER ANGOLI, sistemare ridimensionamento poligoni post taglio
    h0.cutMeshBySegment(mesh,a,b);
    vector<vector<unsigned int>> calculated_vectors; //mi salvo i vettori calcolati
    for (unsigned int i = 0; i < mesh.activatedPolygons.size(); i++){
        unsigned int k = mesh.activatedPolygons[i];
        vector<unsigned int> vi;
        cout << endl;
        for (unsigned int j = 0; j < mesh.VerticesCell2Ds[i].size(); j++){
            vi.push_back(mesh.VerticesCell2Ds[k][j]);
            cout << " " << mesh.VerticesCell2Ds[k][j];
        }
        calculated_vectors.push_back(vi);
    }

    bool flag1 = false; //ora è tempo di controllare se i poligoni sono stati salvati correttamente
    vector<vector<unsigned int>> true_vectors;
    vector<unsigned int> v1 = {4,6,7,0,4,6,7};
    vector<unsigned int> v2 = {4,1,8,6,4,1,8};
    vector<unsigned int> v3 = {6,8,2,5,6,8,2};
    vector<unsigned int> v4 = {7,6,5,3,7,6,5};
    true_vectors.push_back(v1);
    true_vectors.push_back(v2);
    true_vectors.push_back(v3);
    true_vectors.push_back(v4);
    bool different_number_of_polygons = true;
    // if (calculated_vectors.size() != (v1.size()+1)/2){
    //     ASSERT_TRUE(false == different_number_of_polygons);
    // }
    for (unsigned int l = 0; l < 4; l++){ //controllo che i poligoni salvati siano quelli che mi aspetto in senso antiorario
        for (unsigned int i = 0; i < calculated_vectors.size(); i++){
            cout << endl;
            for (unsigned int j = 0; j < calculated_vectors[i].size(); j++){
                //cout << "   " << calculated_vectors[i][j];
                for (unsigned int k = 0; k < (true_vectors[l].size()+1)/2; k++){
                    //cout << "calc: " << calculated_vectors[i][k] << "   v1: " << true_vectors[l][k+j] << endl; //a meno di permutazioni
                    if (calculated_vectors[i][k] != true_vectors[l][k+j]){
                        flag1 = false;
                        break;
                    } else {
                        flag1 = true;
                    }
                }
                if (flag1 == true){
                    break;
                }
            }
            if (flag1 == true){
                break;
            }
        }
        if (flag1 == false){ //se per il primo poligono non ho trovasto corrispondenza nei poligoni calcolati allora ho fallito il test
            break;
        }
    }
    ASSERT_TRUE(flag1 == true); //se è vero allora tutti i poligoni calcolati sono come me li aspetto (giusti e in senso antiorario)
}
//FARE ALTRO ESEMPIO CON TRACCIA CHE TAGLIA I BORDI
