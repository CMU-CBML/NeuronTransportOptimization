#include <iostream>
#include <vector>
#include <array>
#include "BasicDataStructure.h"
#include "NS_3Dsteady.h"
#include "UserSetting.h"
#include "NS_2Dsteady.h"
#include "L2Projection.h"
#include "UserSetting2D.h"
#include "petscmat.h"
#include <sstream>
#include <iomanip>
#include "time.h"
using namespace std;

// static int problem_dim = 3;

static char help[] = "Solve steady Navier Stokes Equation\n";

int main(int argc, char **argv)
{
	stringstream ss_in, ss_out, stmp;
	string path_in, path_out;
	char fld_in[PETSC_MAX_PATH_LEN], fld_out[PETSC_MAX_PATH_LEN];
	PetscBool flg_in, flg_out, flg_dim;
	PetscInt problem_dim;

	int rank, nProcs;
	PetscErrorCode ierr;
	/// start up petsc
	ierr = PetscInitialize(&argc, &argv, (char *)0, help);
	if (ierr)
		return ierr;
	ierr = PetscOptionsGetString(NULL, NULL, "-i", fld_in, PETSC_MAX_PATH_LEN, &flg_in);
	// ierr = PetscOptionsGetString(NULL, NULL, "-o", fld_out, PETSC_MAX_PATH_LEN, &flg_out);
	// if (!flg_in || !flg_out)
	// 	SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate binary file with the -i -o option");
	if (!flg_in)
		SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate working directory with the -i option");
	ierr = PetscOptionsGetInt(NULL, NULL, "-d", &problem_dim, &flg_dim);
	// ierr = PetscOptionsGetString(NULL, NULL, "-o", fld_out, PETSC_MAX_PATH_LEN, &flg_out);
	// if (!flg_in || !flg_out)
	// 	SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate binary file with the -i -o option");
	if (!flg_dim)
		SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate problem dimension with the -d option");
	ss_in << fld_in;
	ss_in >> path_in;
	ss_out << fld_out;
	ss_out >> path_out;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

	// ! NS 3D steady Problem
	if (problem_dim == 3)
	{
		int dimension(3), n_bzmesh;
		vector<double> var;
		vector<array<double, 3>> velocity_node;
		vector<Vertex3D> cpts;
		vector<Element3D> tmesh;
		vector<vector<int>> ele_process;
		vector<double> Vel0, Pre0;
		time_t t0, t1;
		ele_process.resize(nProcs);

		/// Set simulation parameters and mesh
		string fn_mesh(path_in + "controlmesh_label.vtk");
		string fn_bz(path_in + "bzmeshinfo.txt.epart." + to_string(nProcs));
		string fn_velocity(path_in + "initial_velocityfield.txt");
		string fn_parameter(path_in + "simulation_parameter_NS.txt");

		UserSetting *user = new UserSetting;
		user->SetVariables(fn_parameter, var);
		user->ReadMesh(fn_mesh, cpts, tmesh);
		user->ReadVelocityField(fn_velocity, cpts.size(), velocity_node);
		user->AssignProcessor(fn_bz, n_bzmesh, ele_process);
		user->SetInitialCondition(cpts.size(), Vel0, Pre0, cpts, velocity_node, dimension);

		///NS 3D steady Problem
		// Solve Problem
		NS_3Dsteady *NavierStokes3d = new NS_3Dsteady;
		NavierStokes3d->InitializeProblem(cpts.size(), n_bzmesh, Vel0, Pre0, var);
		NavierStokes3d->AssignProcessor(ele_process);
		NavierStokes3d->Run(cpts, tmesh, velocity_node, path_in);

		delete user;
		delete NavierStokes3d;
	}
	else if (problem_dim == 2) // ! NS 2D steady Problem
	{
		int dimension(2), n_bzmesh;
		vector<double> var;
		vector<array<double, 2>> velocity_node;
		vector<Vertex2D> cpts;
		vector<Element2D> tmesh;
		vector<vector<int>> ele_process;
		vector<double> Vel0, Pre0;
		time_t t0, t1;
		ele_process.resize(nProcs);

		/// Set simulation parameters and mesh
		string fn_mesh(path_in + "controlmesh_label.vtk");
		string fn_bz(path_in + "bzmeshinfo.txt.epart." + to_string(nProcs));
		string fn_velocity(path_in + "initial_velocityfield.txt");
		string fn_parameter(path_in + "simulation_parameter_NS.txt");

		UserSetting2D *user = new UserSetting2D;
		user->SetVariables(fn_parameter, var);
		user->ReadMesh(fn_mesh, cpts, tmesh);
		user->ReadVelocityField(fn_velocity, cpts.size(), velocity_node);
		user->AssignProcessor(fn_bz, n_bzmesh, ele_process);
		user->SetInitialCondition(cpts.size(), Vel0, Pre0, cpts, velocity_node, dimension);

		// L2Projection *L2Proj2d = new L2Projection;
		// L2Proj2d->AssignProcessor(ele_process);
		// L2Proj2d->RunL2Projection(cpts, tmesh, path_in);

		//NS 2D steady Problem
		// Solve Problem
		NS_2Dsteady *NavierStokes2d = new NS_2Dsteady;
		NavierStokes2d->InitializeProblem(cpts.size(), n_bzmesh, Vel0, Pre0, var);
		NavierStokes2d->AssignProcessor(ele_process);
		NavierStokes2d->Run(cpts, tmesh, velocity_node, path_in);

		delete user;
		// delete L2Proj2d;
		delete NavierStokes2d;
	}

	PetscPrintf(PETSC_COMM_WORLD, "Done!\n");
	ierr = PetscFinalize();
	CHKERRQ(ierr);

	return 0;
}
