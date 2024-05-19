#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Eigen/Eigen"
#include "classes/Fracture.hpp"
#include "Utils.hpp"

using namespace std;

namespace Utils{

vector<Eigen::Vector3d> calculateDistinctPoints(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b) {
    // la funzione assume che  a e b siano due coppie di punti distinti tra di loro
    vector<Eigen::Vector3d> result;
    result.push_back(a[0]);
    for (const Eigen::Vector3d& b_el : b) {
        if ((a[0]-b_el).norm() >= 8*numeric_limits<double>::epsilon()) {
            result.push_back(b_el);
        }
    }
    if ((a[1]-b[0]).norm() >= 8*numeric_limits<double>::epsilon() && (a[1]-b[1]).norm() >= 8*numeric_limits<double>::epsilon()) {
        result.push_back(a[1]);
    }
    return result;
}

bool compareSegments(vector<Eigen::Vector3d>& a, vector<Eigen::Vector3d>& b) {
    if (a[0] == b[0] && a[1] == b[1]) {
        return true;
    }
    if (a[0] == b[1] && a[1] == b[0]) {
        return true;
    }
    return false;
}

vector<Fracture> fractureInput(const string& filename, array<double, 6>& domain_borders){
    vector<Fracture> output;
    ifstream infile(filename);
    if (infile.fail()) {
        cerr << "Errore nell'apertura del file" << endl;
        return output;
    }

    string line;
    getline(infile,line);
    getline(infile,line);
    unsigned int n_fractures = stoi(line);

    output.reserve(n_fractures);

    double x_coord_min = numeric_limits<double>::max();
    double x_coord_max = numeric_limits<double>::min();
    double y_coord_min = numeric_limits<double>::max();
    double y_coord_max = numeric_limits<double>::min();
    double z_coord_min = numeric_limits<double>::max();
    double z_coord_max = numeric_limits<double>::min();

    for (unsigned int i= 0; i < n_fractures; i++){
        getline(infile,line);

        getline(infile,line);
        istringstream s(line);
        string data;
        getline(s,data,';');
        unsigned int id = stoi(data);
        getline(s,data);
        unsigned int n_vertices = stoi(data);

        getline(infile, data);
        Eigen::MatrixXd vertices;
        vertices.resize(3, n_vertices);

        getline(infile,line);
        istringstream x_line(line);
        for (unsigned int j=0; j<n_vertices; j++) {
            string x_coord;
            getline(x_line, x_coord, ';');
            vertices(0, j) = stod(x_coord);
            if (stod(x_coord) < x_coord_min) {
                x_coord_min = stod(x_coord);
            }
            if (stod(x_coord) > x_coord_max) {
                x_coord_max = stod(x_coord);
            }
        }

        getline(infile,line);
        istringstream y_line(line);
        for (unsigned int j=0; j<n_vertices; j++) {
            string y_coord;
            getline(y_line, y_coord, ';');
            vertices(1, j) = stod(y_coord);
            if (stod(y_coord) < y_coord_min) {
                y_coord_min = stod(y_coord);
            }
            if (stod(y_coord) > y_coord_max) {
                y_coord_max = stod(y_coord);
            }
        }

        getline(infile,line);
        istringstream z_line(line);
        for (unsigned int j=0; j<n_vertices; j++) {
            string z_coord;
            getline(z_line, z_coord, ';');
            vertices(2, j) = stod(z_coord);
            if (stod(z_coord) < z_coord_min) {
                z_coord_min = stod(z_coord);
            }
            if (stod(z_coord) > z_coord_max) {
                z_coord_max = stod(z_coord);
            }
        }

        Fracture element = Fracture(id, n_vertices, vertices);
        output.push_back(element);
    }

    domain_borders[0] = x_coord_min;
    domain_borders[1] = x_coord_max;
    domain_borders[2] = y_coord_min;
    domain_borders[3] = y_coord_max;
    domain_borders[4] = z_coord_min;
    domain_borders[5] = z_coord_max;

    infile.close();
    return output;
}

}
