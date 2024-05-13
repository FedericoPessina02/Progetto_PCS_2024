#include <iostream>
#include <vector>
#include "Utils.hpp"
#

using namespace std;

int main()
{
    vector<Fracture> fractures = Utils::fractureInput("C:/Progetto_PCS_2024/Project/DFN/FR3_data.txt"); //elenco fratture
    TracesMesh meshprova;
    Fracture pippo = fractures[0];
    Fracture other = fractures[1];
    // for (unsigned int i = 0; i < 3; i++){
    //     cout << pippo.normal[i] << "  " << other.normal[i] << endl;
    // }
    // pippo.generateTrace(other,meshprova);

    return 0;
}
