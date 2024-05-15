#include <iostream>
#include <vector>
#include "Utils.hpp"

using namespace std;

int main()
{
    vector<Fracture> fractures = Utils::fractureInput("DFN/FR3_data.txt");
    return 0;
}
