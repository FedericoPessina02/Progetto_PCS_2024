#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include "Utils.hpp"
#include "Algorithms.hpp"
#include "Tests.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    int partition_dimension = 2;
    array<double, 6> domain_borders;
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR3_data.txt", domain_borders); //elenco fratture
    map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, partition_dimension);
    TracesMesh mesh;
    Algorithms::cutTraces(id_to_fractures, mesh);
    Utils::Stampa1("results1.csv",mesh);
    Algorithms::ordinaFract(id_to_fractures, mesh,"results2.csv");
    vector<PolygonalMesh> polygons = Algorithms::cutPolygonalMesh(id_to_fractures, mesh);
    //::testing::InitGoogleTest(&argc,argv);
    // return RUN_ALL_TESTS();
    return 0;
}
