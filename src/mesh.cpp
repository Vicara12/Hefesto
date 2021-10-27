#include "mesh.h"

#include <iostream>  // delete after debug -------------------------->>>>>>>>>>>
using namespace std;

Mesh::Mesh (const tMeshData *mesh) :
        n_volumes(mesh->n_volms), n_boundaries(mesh->n_boundaries)
{
    node = new Volume*[n_volumes+n_boundaries];

    // build mesh of solid volumes
    for (int i = 0; i < n_volumes; i++)
    {
        double volume = mesh->volms_data[i][0];
        double lambda = mesh->volms_data[i][1];
        double qv = mesh->volms_data[i][2];
        const double *surfaces = &(mesh->volms_data[i][3]);
        const Volume *boundaries [PROBLEM_DIM*2];
        double position [PROBLEM_DIM];
        
        for (int j = 0; j < PROBLEM_DIM; j++)
            position[j] = mesh->pos_nodes[i][j];

        // for each boundary index, get a pointer to it
        for (int j = 0; j < PROBLEM_DIM*2; j++)
        {
            int boundary_index = int(mesh->volms_data[i][3+2*PROBLEM_DIM+j]);
            boundaries[j] = node[boundary_index];
        }
            
        node[i] = new SolidVolume(volume, lambda, qv, surfaces, i, position);
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
        const Volume *boundaries [PROBLEM_DIM*2];

        // for each boundary index, get a pointer to it
        for (int j = 0; j < PROBLEM_DIM*2; j++)
        {
            int boundary_index = int(mesh->volms_data[i][3+2*PROBLEM_DIM+j]);
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


void Mesh::solveMesh (DoubleVector &T)
{
    DoubleMatrix eq_sys(n_volumes, DoubleVector(n_volumes+1, 0));

    for (int i = 0; i < n_volumes; i++)
        ((SolidVolume*)node[i])->getEquation(eq_sys[i]);
    
    for (auto l : eq_sys)
    {
        for (auto i : l)
            cout << i << "   ";
        
        cout << endl;
    }
}


void Mesh::solveTransitory (const DoubleVector &T0, DoubleMatrix &T,
                            int time_steps, double t, int store_each)
{
    //
}


void Mesh::printMesh (int from, int to, bool only_volumes)
{
    if (to < 0)
        to = n_volumes + (only_volumes ? 0 : n_boundaries);
    
    for (int i = from; i < to; i++)
        node[i]->print(i);
}


Mesh::~Mesh ()
{
    delete [] node;
}