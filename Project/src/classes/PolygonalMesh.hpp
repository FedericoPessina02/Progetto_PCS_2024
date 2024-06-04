#pragma once
#include "Eigen/Eigen"
#include "vector"

using namespace std;

struct PolygonalMesh {
    unsigned int FractureId;
    Eigen::Vector3d normal;
    vector<unsigned int> IdCell0Ds; // id punti
    vector<Eigen::Vector3d> CoordinateCell0Ds; // coordinate punti
    vector<unsigned int> IdCell1Ds; // id lati
    vector<array<unsigned int,2>> VerticesCell1Ds; // riferimenti agli id dei vertici dei lati
    vector<unsigned int> IdCell2Ds; // id poligoni
    vector<vector<unsigned int>> VerticesCell2Ds; // riferimenti agli id dei vertici dei poligoni
    vector<vector<unsigned int>> EdgesCell2Ds;  // riferimenti agli id dei lati dei poligoni
    vector<unsigned int> activatedPolygons; // vettore di poligoni in uso
    PolygonalMesh() = default;
    unsigned int addPoint(Eigen::Vector3d point) {
        for (unsigned int i = 0; i < CoordinateCell0Ds.size(); i++) {
            if ((point-CoordinateCell0Ds[i]).squaredNorm() < 5*numeric_limits<double>::epsilon()) {
                return i;
            }
        }
        unsigned int new_id = IdCell0Ds.size();
        IdCell0Ds.push_back(new_id);
        CoordinateCell0Ds.push_back(point);
        return new_id;
    }
    unsigned int addEdge(unsigned int a, unsigned int b) {
        for (unsigned int i = 0; i < VerticesCell1Ds.size(); i++) {
            if ((CoordinateCell0Ds[VerticesCell1Ds[i][0]] - CoordinateCell0Ds[a]).squaredNorm() < 5*numeric_limits<double>::epsilon() &&
                (CoordinateCell0Ds[VerticesCell1Ds[i][1]] - CoordinateCell0Ds[b]).squaredNorm() < 5*numeric_limits<double>::epsilon()) {
                return i;
            }
            if ((CoordinateCell0Ds[VerticesCell1Ds[i][1]] - CoordinateCell0Ds[a]).squaredNorm() < 5*numeric_limits<double>::epsilon() &&
                (CoordinateCell0Ds[VerticesCell1Ds[i][0]] - CoordinateCell0Ds[b]).squaredNorm() < 5*numeric_limits<double>::epsilon()) {
                return i;
            }
        }
        unsigned int new_id = IdCell1Ds.size();
        IdCell1Ds.push_back(new_id);
        VerticesCell1Ds.push_back(array<unsigned int,2> {a, b});
        return new_id;
    }
};
