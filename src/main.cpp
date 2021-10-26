#include <iostream>
#include "mesh.h"
using namespace std;

int main ()
{
	tMeshData mesh;

	mesh.n_volms = 3;
	mesh.n_boundaries = 8;

	mesh.pos_nodes = DoubleMatrix({{0.5, 0.5},  // (00): volume 0
							   	   {3.0, 0.5},  // (01): volume 1
							   	   {4.5, 0.5},  // (02): volume 3
							   	   {0.0, 0.5},  // (03): t left
							   	   {0.5, 0.0},  // (04): conv 0 down
							   	   {0.5, 1.0},  // (05): conv 0 up
							   	   {3.0, 0.0},  // (06): conv 1 down
							   	   {3.0, 1.0},  // (07): conv 1 up
							   	   {4.5, 0.0},  // (08): conv 2 down
							   	   {4.5, 1.0},  // (09): conv 2 up
							   	   {6.0, 0.5}}); // (10): conv 2 right

	mesh.volms_data = DoubleMatrix({{1.0, 3.0, -1.0, 1.0,1.0,1.0,1.0, 3,1,4,5}, // 0
									{2.0, 5.0, -2.0, 1.0,1.0,2.0,2.0, 0,2,6,7}, // 1
									{3.0, 7.0, -9.0, 1.0,1.0,3.0,3.0, 1,10,8,9}}); // 3

	mesh.boundary_data = DoubleMatrix({{fixed_T_boundary   , 342.0, 0.0 }, // 3
									   {convection_boundary, 303.0, 32.0}, // 3
									   {convection_boundary, 302.0, 33.0}, // 3
									   {convection_boundary, 306.0, 34.0}, // 3
									   {convection_boundary, 308.0, 35.0}, // 3
									   {convection_boundary, 307.0, 36.0}, // 3
									   {convection_boundary, 309.0, 37.0}, // 3
									   {convection_boundary, 300.0, 38.0}}); // 3

	Mesh new_mesh(&mesh);

	cout << "GENERATED MESH" << endl << endl;

	new_mesh.printMesh(0, -1);

	cout << endl << endl;

	DoubleVector T(mesh.n_volms);

	new_mesh.solveMesh(T);
}