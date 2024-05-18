#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"

using namespace std;

namespace Algorithms {

void cutTraces(vector<Fracture>& fractures, TracesMesh& mesh) {
    for (unsigned int i = 0; i < fractures.size() - 1; i ++) {
        for (unsigned int j = i+1; j<fractures.size(); j++) {
            Fracture& a = fractures[i];
            Fracture& b = fractures[j];
            if ((a.barycenter-b.barycenter).squaredNorm() < 2 * (a.radius + b.radius)) {
                a.generateTrace(b, mesh);
            }
        }
    }
}

}
