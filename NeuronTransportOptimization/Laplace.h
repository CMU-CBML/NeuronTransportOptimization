#ifndef DIFFUSION_H
#define DIFFUSION_H
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mkl.h"
#include "mkl_spblas.h"
#include "BasicDataStructure.h"

using namespace std;
using namespace Eigen;

class Laplace
{
public:
	Laplace();
	vector<double> GetCA();
	vector<double> GetNplus();
	vector<double> GetNminus();
	void GaussInfo(int ng);
	void InitializeProblem(vector<Element3D> &bzmesh, vector<array<double, 3>> &velocity_node, const vector<double>& CA0, const vector<double>& N_plus0, const vector<double>& N_minus0, const vector<double>& var, double tstep);
	void ConstructShapeFunction();
	void CoordinateTransfer(double u, double v, double w, double &x, double &y, double &z, const vector<array<double, 3>>& pt);
	void VelocityTransfer(double u, double v, double w, const vector<array<double, 64>> &cmat, const vector<array<double, 3>>& v_node, double v_tmp[3]);
	void SUPGcoefficient(double velocity[3], double s, const double dudx[3][3], double &tau);

	void SUPGcoefficient2(double v[3], double dNdx[bzpt_num][dim], double &tau_supg);
	void Bzmesh2Tmesh_value(Element3D &bzel, const double value_bzpt[bzpt_num], vector<double>& value_cpt);
	void Tmesh2Bzmesh_value(Element3D &bzel, const vector<double>& value_cpt, double value_bzpt[bzpt_num]);
	void Tmesh2Bzmesh_velocity(Element3D &bzel, const vector<array<double, 3>>& v_cpt_node, vector<array<double, 3>>& v_bzpt_node);

	//Using Bezier Points Data Structure
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[bzpt_num], double dNdx[bzpt_num][dim], double dudx[3][3], double& detJ);
	void WeightingFunction(const double velocity[3], const double& tau, const double Nx[bzpt_num], const double dNdx[bzpt_num][dim], double Wx[bzpt_num]);
	void ElementValue(double Nx[bzpt_num], const double value_bz[bzpt_num], double &value);
	void ElementVelocity(double Nx[bzpt_num], const vector<array<double, 3>> v_node, double v_tmp[3]);
	void ElementMassMatrix(double Weight[bzpt_num], double Nx[bzpt_num], double detJ, double EM[bzpt_num][bzpt_num]);
	void ElementConvectionMatrix(double Nx[bzpt_num], double dNdx[bzpt_num][dim], double v[dim], double detJ, double EC[bzpt_num][bzpt_num]);
	void ElementStiffMatrix(double dNdx[bzpt_num][dim], double detJ, double EK[bzpt_num][bzpt_num]);
	void ElementLoadVector(double Fx[bzpt_num], double detJ, double EF[bzpt_num]);
	void Tangent(double EM[bzpt_num][bzpt_num], double EM_plus[bzpt_num][bzpt_num], double EM_minus[bzpt_num][bzpt_num], double EK[bzpt_num][bzpt_num], double EC_plus[bzpt_num][bzpt_num], double EC_minus[bzpt_num][bzpt_num], double EMatrixSolve[bzpt_num * 2][bzpt_num * 2]);
	void Residual(const double CA_bz[bzpt_num], const double Nplus_bz[bzpt_num], const double Nminus_bz[bzpt_num], double EM[bzpt_num][bzpt_num], double EM_plus[bzpt_num][bzpt_num], double EM_minus[bzpt_num][bzpt_num], double EVectorSolve[bzpt_num * 2]);
	void Assembly(vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve, const vector<int>& IEN, vector<Triplet<double>>& trilist_MatrixSolve, vector<double>& VectorSolve);
	void BuildLinearSystem(const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<array<double, 3>> &cpts, const vector<int>& pid_loc, const vector<int>& label, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus, SparseMatrix<double, RowMajor>& MatrixSolve, vector<double>& VectorSolve);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve);
	void Solver(const vector<array<double, 3>>& cpts, const vector<array<double, 3>> velocity_node, const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<int>& label, const vector<int>& pid_loc);

	void Bzmesh2Tmesh_Matrix(const vector<array<double, 64>>& cmat, const vector<int>& IEN, double EMatrixSolve[bzpt_num * 2][bzpt_num * 2], vector<vector<double>> &EMatrixSolve1);
	void Bzmesh2Tmesh_Vector(const vector<array<double, 64>>& cmat, const vector<int>& IEN,double EVectorSolve[bzpt_num * 2], vector<double> &EVectorSolve1);

	//Using Control Points Data Structure
	void SUPGcoefficientIGA(double s, double v[3], vector<array<double, 3>>& dNdx, double &tau_supg);
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, double dudx[3][3], double& detJ);
	void WeightingFunction(const double velocity[3], const double& tau, const vector<double> &Nx, const vector<array<double, 3>> &dNdx, vector<double> &Wx);
	void ElementValue(const vector<double> &Nx, const vector<double> value_node, double &value);
	void ElementVelocity(const vector<double> &Nx, const vector<array<double, 3>>& v_node, double v_tmp[3]);
	void ElementMassMatrix(vector<double> &Weight, vector<double>& Nx, double detJ, vector<vector<double>>& EM);
	void ElementConvectionMatrix(vector<double>& Nx, vector<array<double, 3>>& dNdx, double v[3], double detJ, vector<vector<double>>& EC);
	void ElementStiffMatrix(vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK);	
	void Tangent(const int nen, vector<vector<double>>& EM, vector<vector<double>>& EM_plus, vector<vector<double>>& EM_minus, vector<vector<double>>& EK, vector<vector<double>>& EC_plus, vector<vector<double>>& EC_minus, vector<vector<double>>& EMatrixSolve);
	void Residual(const int nen, const double CA, const double Nplus, const double Nminus, const vector<double> &Nx, const vector<double> &Npx, const vector<double> &Nmx, const double detJ, vector<double> &EVectorSolve);
	void BuildLinearSystemIGA(const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<array<double, 3>> &cpts, const vector<int>& pid_loc, const vector<int>& label, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus, SparseMatrix<double, RowMajor>& MatrixSolve, vector<double>& VectorSolve);
	//void ShapeFunctionFromBezier(double u, double v, double w, double &Nx);
	
	
	void ConcentrationCal(double u, double v, double w, const Element3D& bzel, double& disp, double& detJ);
	void ConcentrationCal_Coupling_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ);
	void VisualizeVTK(const vector<array<double, 3>>& pts, const vector<Element3D>& mesh, string fn);
	void VisualizeVTK_1(const vector<Element3D>& bzmesh, string fn);
	void VisualizeVTK_2(const vector<Element3D>& bzmesh, string fn);
	void PardisoSolver(SparseMatrix<double, RowMajor>& GK, vector<double>& GF, vector<double>& solution);
	void Run(const vector<array<double, 3>>& pts, const vector<array<double, 3>> velocity_node, const vector<Element3D> &tmesh, vector<Element3D>& bzmesh, const vector<int>& label, const vector<int>& pid_loc);
	double test_function(double x, double y, double z);
private:
	int judge;// judge if the matrix have been solved
	ShapeFunction shape;
	vector<double> Gpt;
	vector<double> wght;
	vector<double> CA;
	vector<double> N_plus;
	vector<double> N_minus;
	vector<array<double, 3>> vplus, vminus;
	vector<double> GF;
	vector<double> VectorSolve;
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-
	SparseMatrix<double, RowMajor> MatrixSolve;
public:
	PardisoLU<SparseMatrix<double, RowMajor>> solver_MIX;

	double DA;
	double dt;
	int ndt;
};

#endif