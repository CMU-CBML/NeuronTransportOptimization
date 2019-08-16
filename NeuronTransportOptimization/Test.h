#ifndef TEST_H
#define TEST_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"

using namespace std;
using namespace Eigen;

class TEST
{
public:
	TEST();
	void GaussInfo(int ng);
	void InitializeProblem(const vector<double>& CA0, const vector<double>& var,double CAi0, double CAs0, double tstep,double length);
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<int>& IEN, double Nx[8], double dNdx[8][3], double& error,double t, double& detJ);
	double calculate_sum(double t, double x, double y, double z);
	double test_function(double x, double y, double z, double &n0, double &n_plus);
	void errorL2(const vector<Element3D>& mesh, double &result, double t);
	void VisualizeVTK(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn);
	void VisualizeVTK_ERROR(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn);
	void Run(const vector<array<double, 3>>& pts, const vector<Element3D>& bzmesh, const vector<double>& Bv, vector<double>& CA1, const vector<int>& pid_loc, const double time, const int judge, string fn);

	vector<double> Gpt;
	vector<double> wght;
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-, k'+,k'-
	vector<double> CA;
	vector<double> N_plus;
	vector<double> N_minus;
	vector<double> CA_Numerical;
	vector<double> Nplus_Numerical;
	vector<double> error_point;
	vector<double> Nplus_error;
	double CAi;
	double CAs;
	double dt;
	double L;
};

#endif