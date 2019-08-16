#ifndef REACTION_H
#define REACTION_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"

using namespace std;
using namespace Eigen;

class LBO
{
public:
	LBO();
	void GaussInfo(int ng);
	void InitializeProblem(const vector<double>& N_plus0, const vector<double>& N_minus0, const vector<double>& var, double tstep);
	void UpdateCA(vector<double>& CA0);
	void BasisFunction(double u, double v, const vector<array<double,3>>& pt, double Nx[4], double dNdx[4][3], double& detJ);
	void ElementMassMatrix(double Nx[4], double detJ, double EM[4][4]);
	void ElementConvectionMatrix(double Nx[4], double dNdx[4][3], double v[3], double detJ, double EC[4][4]);
	void ElementStiffMatrix(double dNdx[4][3], double detJ, double EK[4][4]);
	void Assembly(double EM[4][4], double EK[4][4], double EM_CA[4][4], const vector<int>& IEN, SparseMatrix<double>& GM, SparseMatrix<double>& G_plus, SparseMatrix<double>& G_minus);
	void BuildLinearSystem(const vector<Element2D>& mesh, const double Vplus, const double Vminus, SparseMatrix<double>& GM, SparseMatrix<double>& GK, SparseMatrix<double>& GM_CA);
	void Solver(SparseMatrix<double>& GM, SparseMatrix<double>& GK, SparseMatrix<double>& GM_CA, vector<double>& Bv, vector<double>& CA);
	void VisualizeVTK(const vector<array<double, 3>>& pts, const vector<Element2D>& mesh, string fn);
	void Run(const vector<array<double, 3>>& pts, const vector<Element2D>& bzmesh, vector<double>& CA, vector<double>& Bv, string fn);
	

	vector<double> Gpt;
	vector<double> wght;
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-
	vector<double> N_plus;
	vector<double> N_minus;
	vector<double> CA_B;
	double dt;
};

#endif