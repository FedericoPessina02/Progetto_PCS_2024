#include <iostream>
#include <vector>
#include <chrono>
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
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR50_data.txt", domain_borders); //elenco fratture
    map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, partition_dimension);
    TracesMesh mesh;
    std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
    Algorithms::cutTraces(id_to_fractures, mesh);
    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    std::cout << "traces: " << std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_begin).count() << endl;
    Utils::Stampa1("results1.csv",mesh);
    Algorithms::ordinaFract(id_to_fractures, mesh,"results2.csv");
    t_begin = std::chrono::steady_clock::now();
    vector<PolygonalMesh> polygons = Algorithms::cutPolygonalMeshMultithread(id_to_fractures, mesh);
    t_end = std::chrono::steady_clock::now();
    std::cout << "polygons: " << std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_begin).count() << endl;
    Utils::ExportSTL("polygonal_mesh.stl", polygons);
    // ::testing::InitGoogleTest(&argc,argv);
    // return RUN_ALL_TESTS();
    return 0;
}
