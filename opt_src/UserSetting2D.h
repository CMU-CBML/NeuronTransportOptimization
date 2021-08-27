#ifndef USERSETTING2D_H
#define USERSETTING2D_H

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
#include "Utils.h"
#include <cmath>
using namespace std;

// const bool ReadIC = false;
// const bool ReadBC = false;
// const bool ReadDesire = true;
// // const bool ReadDesire = false;

// const bool VisualizeIC = true;
// const bool VisualizeBC = true;
// const bool VisualizeDesire = true;

// const int phy_dim = 3;
// const int dim = 2;
// const int degree = 3;
// const int bzpt_num = 16;

// // * Neuron Model equation
// // const int state_num = 7;
// // const int ctrl_num = 4;
// // const int result_num = 7;
// // * Burger's equation
// // const int state_num = 2;
// // const int ctrl_num = 2;
// // const int result_num = 6;
// // * Diffusion equation and Convection-diffusion equation
// const int state_num = 1;
// const int ctrl_num = 1;
// const int result_num = 3;

// const int time_int = 0; // * 0 - steady state; 1 - trapezoidal; 2 - rectangle

// const int debug_rank = 0; // * For Parallel implementation debug

// const double PI = 4 * atan(1.0);

//problem setting
class UserSetting2D
{
private:
	static const bool ReadIC = false;
	static const bool ReadBC = false;
	static const bool ReadDesire = true;
	// static const bool ReadDesire = false;

	static const bool VisualizeIC = true;
	static const bool VisualizeBC = true;
	static const bool VisualizeDesire = true;

	static const int phy_dim = 3;
	static const int dim = 2;
	static const int degree = 3;
	static const int bzpt_num = 16;

	// * Neuron Model equation
	// const int state_num = 7;
	// const int ctrl_num = 4;
	// const int result_num = 7;
	// * Burger's equation
	// const int state_num = 2;
	// const int ctrl_num = 2;
	// const int result_num = 6;
	// * Diffusion equation and Convection-diffusion equation
	static const int state_num = 1;
	static const int ctrl_num = 1;
	static const int result_num = 3;

	static const int time_int = 0; // * 0 - steady state; 1 - trapezoidal; 2 - rectangle
public:
	int n_bzmesh;
	int n_bcpt;
	int n_bcval;

	vector<vector<int>> ele_process;
	vector<vector<int>> bzmeshinfo;
	vector<Element2D> bzmesh_process;
	DM dm_bzmesh;

	vector<double> var;
	vector<Vertex2D> pts;
	vector<Element2D> mesh;
	vector<double> val_ini[2];
	vector<double> val_bc[2];
	vector<double> val_desire[2];
	vector<int> bc_flag;	   // global node index -> index without bc pts
	vector<int> nonbc_mapping; // index without bc pts -> global node index
	vector<int> bc_mapping;	   // index of bc pts -> global node index

	string output_dir;
	string work_dir;

private:
	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;
	int nTstep, nPoint;

	vector<double> Gpt;
	vector<double> wght;
	const double PI = 4 * atan(1.0);

private:
	void SetVariables(string fn_par);
	void SetDesireState(string fn_in, string fn_out);
	void SetInitialCondition(string fn_in, string fn_out);
	void SetBoundaryCondition(string fn_in, string fn_out);
	void SetBoundaryMapping();

	void TXTWriteIC(string fn_out);
	void TXTWriteBC(string fn_out);
	void TXTWriteDesire(string fn_out);

	void VTKVisualizeIC(string fn_out);
	void VTKVisualizeBC(string fn_out);
	void VTKVisualizeDesire(string fn_out);

	void ReadMesh(string fn);
	void Readbzmeshinfo(string fn);
	void CreateDMPlexMesh(string fn);
	void ReadBezierElementProcess(string fn);
	void ReadVelocityField(string fn);
	void AssignProcessor(string fn);

public:
	UserSetting2D();
	~UserSetting2D();
	bool GetReadDesire() const;
	void InitializeUserSetting(string fn_in, string fn_out);
	void DesireStateFunction(double x, double y, double z, double t, double result[state_num]) const;
};
#endif