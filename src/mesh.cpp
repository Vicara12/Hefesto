#include "mesh.h"
#include <iostream>


Mesh::Mesh (const tMeshData *mesh) :
        n_volumes(mesh->n_volms), n_boundaries(mesh->n_boundaries),
        problem_dim_(mesh->problem_dimensions),
        node(mesh->n_volms+mesh->n_boundaries)
{
    // build mesh of solid volumes
    for (int i = 0; i < n_volumes; i++)
    {
        double volume = mesh->volms_data[i][0];
        double lambda = mesh->volms_data[i][1];
        double qv = mesh->volms_data[i][2];
            
        node[i] = new SolidVolume(problem_dim_, volume, lambda, qv,
                            mesh->surface_volumes[i], i, mesh->pos_volumes[i]);
    }

    // initialize boundary objects
    for (int i = 0; i < n_boundaries; i++)
    {
        VType node_type = VType(int(mesh->boundary_data[i][0]));

        // each type of volume needs to be initialized differently
        if (node_type == convection_boundary)
        {
            double T_ext = mesh->boundary_data[i][1];
            double alpha = mesh->boundary_data[i][2];

            node[n_volumes+i] = new ConvectionBoundary(T_ext, alpha);
        }
        else if (node_type == fixed_T_boundary)
        {
            double T = mesh->boundary_data[i][1];
            double distance = mesh->boundary_data[i][2];

            node[n_volumes+i] = new FixedTBoundary(T, distance);
        }
        else
        {
            throw "unrecognized boundary type at mesh assembly";
        }
    }

    // now that all the nodes are already generated, buil conectivities
    for (int i = 0; i < n_volumes; i++)
    {
        std::vector<const Volume*> boundaries(2*problem_dim_);

        // for each boundary index, get a pointer to it
        for (int j = 0; j < problem_dim_*2; j++)
        {
            int boundary_index = int(mesh->connectivity_volumes[i][j]);
            boundaries[j] = node[boundary_index];
        }
            
        ((SolidVolume*)node[i])->setBoundaries(boundaries);
    }
}


int Mesh::getNumVolumes () const
{
    return n_volumes;
}


int Mesh::getNumBoundaries () const
{
    return n_boundaries;
}


void Mesh::setNodeData (const DoubleMatrix &node_data)
{
    //
}


void Mesh::setBoundaryData (const DoubleMatrix &boundary_data)
{
    //
}


double Mesh::solveMesh (void(*solver)(const DoubleMatrix&, DoubleVector&, double, bool),
                      DoubleVector &T, double tolerance, bool check_solution,
                      bool verbose)
{
    DoubleMatrix eq_sys(n_volumes, DoubleVector(n_volumes+1, 0));

    for (int i = 0; i < n_volumes; i++)
        ((SolidVolume*)node[i])->getEquation(eq_sys[i]);
    
    solver(eq_sys, T, tolerance, verbose);

    double max_error = 0;

    if (check_solution)
    {
        for (int i = 0; i < n_volumes; i++)
        {
            double this_error = -eq_sys[i][n_volumes];

            for (int j = 0; j < n_volumes; j++)
                this_error += eq_sys[i][j]*T[j];
            
            if (i == 0 or this_error > max_error)
                max_error = this_error;
        }
    }

    return max_error;
}


void Mesh::solveTransitory (const DoubleVector &T0, DoubleMatrix &T,
                            int time_steps, double t, int store_each)
{
    //
}


void Mesh::printMesh (int from, int to, bool only_volumes) const
{
    if (to < 0)
        to = n_volumes + (only_volumes ? 0 : n_boundaries);
    
    for (int i = from; i < to; i++)
        node[i]->print(i);
}


void Mesh::printNode (int index) const
{
    node[index]->print(index);
}


double Mesh::checkEnergyBalance (const DoubleVector &T) const
{
    double max_value = ((SolidVolume*)node[0])->checkEnergyBalance(T);

    for (int i = 1; i < n_volumes; i++)
    {
        double value = ((SolidVolume*)node[i])->checkEnergyBalance(T);

        if (value > max_value)
            max_value = value;
    }

    return max_value;
}


Mesh::~Mesh ()
{
    for (int i = 0; i < n_boundaries+n_volumes; i++)
        delete node[i];
}