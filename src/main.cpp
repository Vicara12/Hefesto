#include <iostream>
#include <math.h>
#include "mesh.h"
#include "solver.h"
using namespace std;


#define N_ELMS 100
#define SOLVER_TOLERANCE 1e-6


void buildCylindricalFinMesh (tMeshData &mesh)
{
	int n_elms = N_ELMS;
	double Ra = 0.05;
    double Rb = 0.13;
    double e = 0.003;
    double Ta = 200;
    double Tg = 25;
    double lambda = 22;
    double alpha = 100;

	double delta_r = (Rb-Ra)/n_elms;

	mesh.n_volms = n_elms;
	mesh.n_boundaries = 3; // one fixed T for the join with the tube, other for upper
						   // and lower surfaces and other for the adiabatic tip

	mesh.pos_volumes = DoubleMatrix(n_elms, DoubleVector(PROBLEM_DIM, e/2));
	mesh.volms_data = DoubleMatrix(n_elms, DoubleVector(PROBLEM_DIM*4+3, e/2));

	// initialize all data for the volumes
	for (int i = 0; i < n_elms; i++)
	{
		double r = Ra + delta_r*(i + 0.5);
		double vert_surface = M_PI*(pow(r+delta_r/2,2) - pow(r-delta_r/2,2));
		mesh.pos_volumes[i][0] = r;
		mesh.volms_data[i][0] = e*vert_surface;
		mesh.volms_data[i][1] = lambda;
		mesh.volms_data[i][2] = 0;
		mesh.volms_data[i][3] = 2*M_PI*(r-delta_r/2)*e; // left surface
		mesh.volms_data[i][4] = 2*M_PI*(r+delta_r/2)*e; // right surface
		mesh.volms_data[i][5] = vert_surface; // down surface
		mesh.volms_data[i][6] = vert_surface; // up surface

		// joint with the tube
		if (i == 0)
		{
			mesh.volms_data[i][7] = n_elms; // fiex t boundary node
			mesh.volms_data[i][8] = i+1;	// second node
		}
		// tip of the fin
		else if (i == n_elms-1)
		{
			mesh.volms_data[i][7] = i-1; // the volume before the last
			mesh.volms_data[i][8] = n_elms+2; // adiabatic end node
		}
		// intermediate nodes
		else
		{
			mesh.volms_data[i][7] = i-1;	// left volume
			mesh.volms_data[i][8] = i+1;	// right volume
		}

		// upper and lower part are convective boundary
		mesh.volms_data[i][9]  = n_elms+1;
		mesh.volms_data[i][10] = n_elms+1;
	}
	
	// initialize all data for the boundaries
	mesh.boundary_data = DoubleMatrix({{fixed_T_boundary, Ta, delta_r/2}, // n_elms:   fixed T boundary
									   {convection_boundary, Tg, alpha},  // n_elms+1: air convection
									   {convection_boundary, Tg, 0.0}});  // n_elms+2: adiabatic tip

}

void buildTestMesh (tMeshData &mesh)
{
	mesh.n_volms = 3;
	mesh.n_boundaries = 8;

	mesh.pos_volumes = DoubleMatrix({{0.5, 0.5},  // (00): volume 0
							   	    {3.0, 0.5},  // (01): volume 1
							   	    {4.5, 0.5}}); // (02): volume 3

	mesh.volms_data = DoubleMatrix({{1.0, 3.0, -1.0, 1.0,1.0,1.0,1.0, 3,1,4,5}, // 0
									{2.0, 5.0, -2.0, 1.0,1.0,2.0,2.0, 0,2,6,7}, // 1
									{3.0, 7.0, -9.0, 1.0,1.0,3.0,3.0, 1,10,8,9}}); // 3

	mesh.boundary_data = DoubleMatrix({{fixed_T_boundary   , 342.0, 0.5 }, // 3
									   {convection_boundary, 303.0, 32.0}, // 4
									   {convection_boundary, 302.0, 33.0}, // 5
									   {convection_boundary, 306.0, 34.0}, // 6
									   {convection_boundary, 308.0, 35.0}, // 7
									   {convection_boundary, 307.0, 36.0}, // 8
									   {convection_boundary, 309.0, 37.0}, // 9
									   {convection_boundary, 300.0, 0.0}}); // 10
}


int main ()
{
	tMeshData mesh_data;

	buildCylindricalFinMesh(mesh_data);
	//buildTestMesh(mesh_data); 

	Mesh new_mesh(&mesh_data);

	new_mesh.printMesh(0, -1);

	DoubleVector T;

	new_mesh.solveMesh(gaussSeidel, T, 1e-6, false);

	cout << endl << endl;
	cout << "Finished solving mesh, final solution: " << endl << endl;

	for (int i = 0; i < mesh_data.n_volms; i++)
	{
		cout << "Volume " << i << " T = " << T[i] << " K\n";
	}
}