#include <iostream>
#include <vector>
#include <array>
#include "BasicDataStructure.h"
// #include "NS_3Dsteady.h"
// #include "TransportOpt2D.h"
#include "DiffConv2D.h"
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
	PetscErrorCode ierr;
	ierr = PetscInitialize(&argc, &argv, (char *)0, help);
	if (ierr)
		return ierr;

	stringstream ss_in, ss_out, stmp;
	string path_in, path_out;
	char fld_in[PETSC_MAX_PATH_LEN], fld_out[PETSC_MAX_PATH_LEN];
	// ss << argv[1];
	// ss >> path;
	// stmp << argv[2];
	// int n_process = atoi(argv[2]);

	int comRank, comSize;
	int problem_dim;

	PetscBool flg_in, flg_out, flg_dim;
	/// start up petsc option
	ierr = PetscOptionsGetString(NULL, NULL, "-i", fld_in, PETSC_MAX_PATH_LEN, &flg_in);
	ierr = PetscOptionsGetString(NULL, NULL, "-o", fld_out, PETSC_MAX_PATH_LEN, &flg_out);
	ierr = PetscOptionsGetInt(NULL, NULL, "-d", &problem_dim, &flg_dim);
	if (!flg_in || !flg_out)
		SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate I/O path with the -i -o option");
	if (!flg_dim)
		SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate problem dimension with the -d option");
	ss_in << fld_in;
	ss_in >> path_in;
	ss_out << fld_out;
	ss_out >> path_out;

	MPI_Comm_rank(PETSC_COMM_WORLD, &comRank);
	MPI_Comm_size(PETSC_COMM_WORLD, &comSize);

	// printf("Hello world from rank %d out of %d processors\n", comRank, comSize);

	// New UserSetting2D
	if (problem_dim == 2)
	{

		UserSetting2D *ctx = new UserSetting2D;
		DiffConv2D *diffconv2d = new DiffConv2D;

		// TransportOpt2D *transopt2d = new TransportOpt2D;
		ctx->InitializeUserSetting(path_in, path_out);

		// * Diffusion convection test
		diffconv2d->Run(ctx);

		// * Burger's equation test
		// diffconv2d->Run(ctx);

		// * Neuron Transport test
		// transopt2d->Run(ctx);

		PetscPrintf(PETSC_COMM_WORLD, "Done!\n");

		delete ctx;
		delete diffconv2d;
		// delete transopt2d;
	}
	// else if (problem_dim == 3)
	// {
	// 	UserSetting3D *ctx = new UserSetting3D;
	// 	DiffConv3D *diffconv3d = new DiffConv3D;

	// 	ctx->InitializeUserSetting(path_in, path_out);
	// 	diffconv3d->Run(ctx);
	// 	PetscPrintf(PETSC_COMM_WORLD, "Done!\n");

	// 	delete ctx;
	// 	delete diffconv3d;
	// }

	ierr = PetscFinalize();

	return 0;
}
