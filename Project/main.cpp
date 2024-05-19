#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    array<double, 6> domain_borders;
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR362_data.txt", domain_borders); //elenco fratture
    Algorithms::assignPartition(fractures, domain_borders, 2);
    TracesMesh mesh;
    Algorithms::cutTraces(fractures, mesh);
    return 0;
}
