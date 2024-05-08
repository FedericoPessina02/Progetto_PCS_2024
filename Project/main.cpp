#include <iostream>
#include <vector>
#include "Utils.hpp"

using namespace std;

int main()
{
    vector<Fracture> fractures = Utils::fractureInput("C:/Progetto_PCS_2024/Project/DFN/FR3_data.txt");
    return 0;
}
