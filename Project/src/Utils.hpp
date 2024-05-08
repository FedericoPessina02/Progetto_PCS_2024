#pragma once
#include <iostream>
#include <vector>
#include "classes/Fracture.hpp"

using namespace std;

namespace Utils{

vector<Fracture> fractureInput(const string& filename);

};
