#ifndef SOLVER_H_
#define SOLVER_H_

#include "definitions.h"

void gaussSeidel (const DoubleMatrix &system, DoubleVector &solution,
                 double tolerance, bool verbose);

void TDMA (const DoubleMatrix &system, DoubleVector &solution,
           double void_parameter, bool verbose);

#endif