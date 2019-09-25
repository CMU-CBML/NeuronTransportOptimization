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
const bool VisualizeIC = true;
const bool VisualizeBC = true;

//problem setting
class UserSetting2D
{
public:
	int n_bzmesh;
	vector<vector<int>> ele_process;
	vector<double> var;
	vector<Vertex2D> pts;
	vector<Element2D> mesh;
	vector<double> val_ini[18];
	vector<double> val_bc[7];

	string work_dir;


private:
	void SetVariables(string fn_par);
	void SetInitialCondition(string fn_in, string fn_out);
	void SetBoundaryCondition(string fn_in, string fn_out);

	void TXTWriteIC(string fn_out);
	void TXTWriteBC(string fn_out);

	void VTKVisualizeIC(string fn_out);
	void VTKVisualizeBC(string fn_out);
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