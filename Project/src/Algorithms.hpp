#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "vector"
#include "classes/Fracture.hpp"
#include "classes/TracesMesh.hpp"

using namespace std;

namespace Algorithms {

void cutTraces(vector<Fracture> fractures, TracesMesh& mesh);

}
