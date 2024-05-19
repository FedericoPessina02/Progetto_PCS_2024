#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    array<double, 6> domain_borders;
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR200_data.txt", domain_borders); //elenco fratture
    map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, 2);
    TracesMesh mesh;
    Algorithms::cutTraces(id_to_fractures, mesh, 8);
    return 0;
}
