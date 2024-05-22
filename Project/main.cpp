#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    array<double, 6> domain_borders;
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR3_data.txt", domain_borders); //elenco fratture
    map<int, vector<Fracture>> id_to_fractures = Algorithms::assignPartition(fractures, domain_borders, 2);
    TracesMesh mesh;
    Algorithms::cutTraces(id_to_fractures, mesh, 8);
    Utils::Stampa1("results1.csv",mesh);
    Algorithms::ordinaFract(id_to_fractures, mesh,"results2.csv");


    return 0;
}
