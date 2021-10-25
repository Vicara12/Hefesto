#include "mesh.h"

Mesh::Mesh (const tMeshData *mesh)
{
    //
}


int Mesh::getNumNodes () const
{
    //
}


int Mesh::getNumBoundaries () const
{
    //
}


void Mesh::setNodeData (double **node_data)
{
    //
}


void Mesh::setBoundaryData (double **boundary_data)
{
    //
}


void Mesh::solveMesh (double *T)
{
    //
}


void Mesh::solveTransitory (double *T0, double **T, int time_steps,
                            double t, int store_each)
{
    //
}