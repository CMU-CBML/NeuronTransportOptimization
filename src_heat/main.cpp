#include <iostream>
#include <vector>
#include <array>
#include "BasicDataStructure.h"
// #include "NS_3Dsteady.h"
#include "DiffConv2D.h"
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

	stringstream ss, stmp;
	string path;
	char file[PETSC_MAX_PATH_LEN];
	// ss << argv[1];
	// ss >> path;
	// stmp << argv[2];
	// int n_process = atoi(argv[2]);

	int rank, nProcs;
	PetscErrorCode ierr;
	PetscBool flg;

	/// start up petsc
	ierr = PetscInitialize(&argc, &argv, (char *)0, help);
	if (ierr)
		return ierr;

	ierr = PetscOptionsGetString(NULL, NULL, "-f", file, PETSC_MAX_PATH_LEN, &flg);
	if (!flg)
		SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate binary file with the -f option");

	ss << file;
	ss >> path;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

	/// New UserSetting2D
	UserSetting2D *ctx = new UserSetting2D;
	DiffConv2D *diffconv2d = new DiffConv2D;
	TransportOpt2D *transopt2d = new TransportOpt2D;
	ctx->InitializeUserSetting(path);

	// Diffusion convection test
	diffconv2d->Run(ctx);

	// Burger's equation test
	// diffconv2d->Run(ctx);

	// Neuron Transport test
	// transopt2d->Run(ctx);

	PetscPrintf(PETSC_COMM_WORLD, "Done!\n");

	delete ctx;
	delete diffconv2d;
	delete transopt2d;

	ierr = PetscFinalize();
	CHKERRQ(ierr);

	return 0;
}
