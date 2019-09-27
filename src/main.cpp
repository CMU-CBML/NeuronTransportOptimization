#include <iostream>
#include <vector>
#include <array>
#include "BasicDataStructure.h"
// #include "NS_3Dsteady.h"
#include "TransportOpt2D.h"
#include "UserSetting2D.h"
#include <sstream>
#include <iomanip>
#include "time.h"
#include <petscsys.h>
#include <petsc.h>



using namespace std;

static char help[] = "Solve 2D Neuron Transport Opt Problem\n";

int main(int argc, char **argv)
{

	if (argc == 3) 
	{
		stringstream ss, stmp;
		string path;
		ss << argv[1];
		ss >> path;
		stmp << argv[2];
		int n_process = atoi(argv[2]);

		int rank, nProcs;
		PetscErrorCode ierr;
		/// start up petsc
		ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

		int dim(2), n_bzmesh;
    	int ntstep(3);

		vector<double> var;		
		vector<array<double, 2>>velocity_node;
		vector<Vertex2D> cpts;
		vector<Element2D> tmesh;
		vector<vector<int>> ele_process;
		vector<double> Vel0, Pre0;
		time_t t0, t1;
		ele_process.resize(nProcs);

		/// New UserSetting2D
		UserSetting2D *ctx = new UserSetting2D;
		TransportOpt2D *transopt2d = new TransportOpt2D;
		
		ctx -> InitializeUserSetting(path);
		transopt2d -> Run(ctx);

		//transopt2d -> InitializeProblem(ctx);

		// transopt2d->InitializeProblem(cpts.size(), n_bzmesh, Vel0, Pre0, var);
		// transopt2d->AssignProcessor(ele_process);
		//transopt2d->Run(cpts, tmesh, velocity_node, path);
		
		// /// Original code
		// { /// Set simulation parameters and mesh
		// 	string fn_mesh(path + "controlmesh.vtk");
		// 	string fn_bz(path + "bzmeshinfo.txt.epart." + stmp.str());
		// 	string fn_velocity(path + "initial_velocityfield.txt");
		// 	string fn_parameter(path + "simulation_parameter.txt");

		// 	UserSetting2D *user = new UserSetting2D;
		// 	user->SetVariables(fn_parameter, var);
		// 	user->ReadMesh(fn_mesh, cpts, tmesh);
		// 	user->ReadVelocityField(fn_velocity, cpts.size(), velocity_node);
		// 	user->AssignProcessor(fn_bz, n_bzmesh, ele_process);
		// 	user->SetInitialCondition(cpts.size(), Vel0, Pre0, cpts, velocity_node);

		// 	///NS 3D steady Problem
		// 	// Solve Problem

		// 	cout << "Setup Simulation" << endl;
		// 	cout << "\n \n";

		// 	TransportOpt2D *transopt2d = new TransportOpt2D;
		// 	transopt2d->InitializeProblem(cpts.size(), n_bzmesh, Vel0, Pre0, var);
		// 	transopt2d->AssignProcessor(ele_process);
		// 	transopt2d->Run(cpts, tmesh, velocity_node, path);
		// }
		PetscPrintf(PETSC_COMM_WORLD, "Done!\n");

		delete ctx;
		delete transopt2d;

		ierr = PetscFinalize(); CHKERRQ(ierr);
	}
	else if (argc > 3) {
		cout << "Too many arguments.\n";
	}
	else {
		cout << "Two argument expected.\n";
	}
	return 0;
}
