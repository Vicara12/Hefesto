#ifndef MESH_H_
#define MESH_H_

#include "volume.h"


typedef struct _tMeshData
{
    int problem_dimensions;
    int n_nodes;
    int n_boundaries;

    // 2D array with the position of each volume in each row
    double **pos_nodes;
    double **pos_boundaries; // position for convection boundaries is not needed
                             // but must be included (the value given doesn't matter)

    // for each node: volume, lambda, qv, boundary surfaces value (2*problem dimension values),
    //       boundary volumes (2*problem dimension values), position (2*problem dimension values)
    double **node_data;

    // for each boundary: type (VType value), T_ext/T, alpha/0
    double **boundary_data;
} tMeshData;


class Mesh
{
public:

    Mesh (const tMeshData *mesh);

    int getNumNodes () const;
    int getNumBoundaries () const;

    // node_data and boundary_data must follow the same format as the one
    // from tMeshData and be of corresponding lenght
    void setNodeData (double **node_data);
    void setBoundaryData (double **boundary_data);

    // T must be a 1D array of length n_nodes. The ith position contains the
    // temperature of the ith node.
    void solveMesh (double *T);

    // T0 must be a 1D array of length n_nodes detailing the initial conditions,
    // T must be a 2D array of size time_steps/store_each x n_nodes and it stores
    // the temperature of all the nodes each store_each timesteps, time_steps is
    // the quantity of time iterations from 0 to t, t is the total simulation
    // time and store_each indicates each how many time_steps the temperature
    // data is stored
    void solveTransitory (double *T0,
                          double **T,
                          int time_steps,
                          double t,
                          int store_each);

private:

    //
};

#endif