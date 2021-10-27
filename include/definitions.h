#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <vector>

#define PROBLEM_DIM 2

enum VType {solid, convection_boundary, fixed_T_boundary};

typedef std::vector<std::vector<double>> DoubleMatrix;
typedef std::vector<double> DoubleVector;

#endif