#include <iostream>
#include "solver.h"


void gaussSeidel (const DoubleMatrix &system, DoubleVector &solution,
                 double tolerance, bool verbose)
{
    int n_nodes = system.size();
    solution = DoubleVector(n_nodes, 0);
    double max_error = tolerance+1;

    if (verbose)
        std::cout << "Beggining Gauss-Seidel" << std::endl;
    
    int n_iter = 0;

    while (max_error > tolerance)
    {
        max_error = 0;

        for (int i = 0; i < n_nodes; i++)
        {
            double value = system[i][n_nodes];

            for (int j = 0; j < n_nodes; j++)
                if (i != j)
                    value -= system[i][j]*solution[j];
            
            value = value/system[i][i];

            double current_error = value - solution[i];

            if (current_error < 0)
                current_error *= -1;

            if (current_error > max_error)
                max_error = current_error;
                        
            solution[i] = value;
        }

        if (verbose)
            std::cout << " - Iteration " << n_iter <<" error: "
                      << max_error << std::endl;
    
        n_iter++;
    }
}


void TDMA (const DoubleMatrix &system, DoubleVector &solution,
           double void_parameter, bool verbose)
{
    throw "NOT IMPLEMENTED";
}