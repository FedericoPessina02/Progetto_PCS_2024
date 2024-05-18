#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"

using namespace std;

namespace Algorithms {

void cutTraces(vector<Fracture>& fractures, TracesMesh& mesh) {
    fractures[0].generateTrace(fractures[1], mesh);
    fractures[0].generateTrace(fractures[2], mesh);
    fractures[1].generateTrace(fractures[2], mesh);
}

}
