#include "Fracture.hpp"
#include <iostream>
#include "Eigen/Eigen"
#include "vector"

array<double,3> crossProduct(const array<double,3>& lato1,const array<double,3>& lato2){
    array<double,3> normal;
    normal[0] = lato1[1]*lato2[2]-lato1[2]*lato2[1];
    normal[1] = lato1[2]*lato2[0]-lato1[0]*lato2[2];
    normal[2] = lato1[0]*lato2[1]-lato1[1]*lato2[0];
    return normal;
}
