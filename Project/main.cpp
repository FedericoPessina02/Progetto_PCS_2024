#include <iostream>
#include <vector>
#include "Utils.hpp"
#include "Algorithms.hpp"

using namespace std;

int main()
{
    vector<Fracture> fractures = Utils::fractureInput("C:/PCS_esercitazioni/Progetto_PCS_2024/Project/DFN/FR362_data.txt"); //elenco fratture
    TracesMesh mesh;
    Algorithms::cutTraces(fractures, mesh);
    return 0;
}
