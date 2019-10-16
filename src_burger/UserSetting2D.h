#ifndef USERSETTING2D_H
#define USERSETTING2D_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
#include <cmath>
using namespace std;

const int dim = 2;

const bool ReadIC = false;
const bool ReadBC = false;
const bool ReadDesire = false;

const bool VisualizeIC = true;
const bool VisualizeBC = true;
const bool VisualizeDesire = true;

//problem setting
class UserSetting2D
{
public:
	int n_bzmesh;
	int n_bcpt;
	int n_bcval;

	vector<vector<int>> ele_process;
	vector<double> var;
	vector<Vertex2D> pts;
	vector<Element2D> mesh;
	vector<double> val_ini[2];
	vector<double> val_bc[2];
	vector<double> val_desire[2];
	vector<int> bc_flag; // global node index -> index without bc pts 
	vector<int> nonbc_mapping; // index without bc pts -> global node index
	vector<int> bc_mapping; // index of bc pts -> global node index

	string work_dir;


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
	void ReadVelocityField(string fn);
	void AssignProcessor(string fn);

public:
	UserSetting2D();
	void InitializeUserSetting(string fn);



	void SetVariables(string fn_par, vector<double> &var);
	void SetInitialCondition(int ndof, vector<double> &Vel0, vector<double> &Pre0, vector<Vertex2D> &pts, const vector<array<double, 2>> velocity_node);
	void ReadMesh(string fn, vector<Vertex2D> &pts, vector<Element2D> &mesh);
	void ReadVelocityField(string fn, int npts, vector<array<double, 2>> &velocity);
	void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);
};
#endif