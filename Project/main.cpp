#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    int partition_dimension = 2;
    array<double, 6> domain_borders;
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR3_data.txt", domain_borders); //elenco fratture
    map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, partition_dimension);
    TracesMesh mesh;
    cout << id_to_fractures.size();
    Algorithms::cutTraces(id_to_fractures, mesh);
    vector<PolygonalMesh> polygons = Algorithms::cutPolygonalMesh(id_to_fractures);
    return 0;
}
