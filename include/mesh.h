#ifndef MESH_H_
#define MESH_H_

#include "volume.h"
#include "definitions.h"


typedef struct _tMeshData
{
    int n_volms;
    int n_boundaries;

    // 2D array with the position of each volume and boundary. First all volumes
    // must be listed and then the boudnaries. The id of a node is denoted
    // by its position in this array
    // position for convection boundaries is not needed
    // but must be included (the value given doesn't matter)
    DoubleMatrix pos_volumes;

    // for each node: volume, lambda, qv, boundary surfaces value (2*problem dimension values),
    //       neighbour nodes (2*problem dimension values)
    DoubleMatrix volms_data;

    // for each boundary: type (VType value), T_ext/T, alpha/distance
    // note: volume types (VType) are defined in definitions.h
    DoubleMatrix boundary_data;
} tMeshData;


class Mesh
{
public:

    Mesh (const tMeshData *mesh);

    int getNumVolumes () const;
    int getNumBoundaries () const;

    // node_data and boundary_data must follow the same format as the one
    // from tMeshData and be of corresponding lenght
    void setNodeData (const DoubleMatrix &node_data);
    void setBoundaryData (const DoubleMatrix &boundary_data);

    // The ith position of T contains the temperature of the ith volume.
    // T must be a 1D vector of size n_nodes and solver the name of a solver from
    // solver.h. If the selected solver has machine precission and no tolerance
    // is needed, a dummy value must be provided anyway.
    //  If verbose = true, the solver will output information about the progress.
    void solveMesh (void(*solver)(const DoubleMatrix&, DoubleVector&, double, bool),
                    DoubleVector &T, double tolerance, bool verbose = false);

    // T0 must be a 1D vector of length n_nodes detailing the initial conditions,
    // T must be a 2D vector of size time_steps/store_each x n_nodes and
    // it stores the temperature of all the nodes each store_each timesteps,
    // time_steps is the quantity of time iterations from 0 to t, t is the total
    // simulation time and store_each indicates each how many time_steps the
    // temperature data is stored
    void solveTransitory (const DoubleVector &T0,
                          DoubleMatrix &T,
                          int time_steps,
                          double t,
                          int store_each);
    
    // from and to indicate the node index that will be printed
    // set to = -1 to print until the last node
    void printMesh (int from, int to, bool only_volumes = true) const;

    void printNode (int index) const;

    // Checks if the solution satisfies energy balance for the current mesh
    double checkEnergyBalance (const DoubleVector &T);

    ~Mesh ();
private:

    int n_volumes;
    int n_boundaries;

    Volume **node; // array with all the nodes
};

#endif