#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    vector<Fracture> fractures = Utils::fractureInput("./DFN/FR10_data.txt"); //elenco fratture
    TracesMesh mesh;
    Algorithms::cutTraces(fractures, mesh);
    return 0;
}
