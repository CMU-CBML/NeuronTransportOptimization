#ifndef NS_2DSTEADY_H
#define NS_2DSTEADY_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include "BasicDataStructure.h"
#include "Utils.h"
#include "time.h"

using namespace std;

// const int degree = 3;
// const int dim = 2;
// const int bzpt_num = 16;
// const int dof_all = 3;

class NS_2Dsteady
{
private:
	static const int degree = 3;
	static const int dim = 2;
	static const int phy_dim = 3;
	static const int bzpt_num = 16;
	static const int dof_all = 3;

	vector<double> Gpt;
	vector<double> wght;
	const double PI = 4.0 * atan(1.0);

	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;

	int rstart, rend;
	int n_bzmesh;
	vector<int> ele_process;
	vector<Element2D> bzmesh_process;

	KSP ksp;
	PC pc;
	Mat GK;
	Vec GR;
	Vec temp_solution;

	double dt;
	double velocity_max;
	double nu;
	double rou;
	double alphaM;
	double alphaF;
	double Gama;
	double fx, fy, fz; //elemental body force

	vector<double> Vel;
	vector<double> Pre;
	vector<double> gh;	//prescribed boundary velocities
	vector<double> par; //Dn0, v_plus, v_minus, k+, k-,k'+,k'-

public:
	NS_2Dsteady();

private:
	void BodyForce(double x, double y, double z, double &Fx, double &Fy, double &Fz);
	void ReadBezierElementProcess(string fn);

	/*Analysis*/
	void GaussInfo(int ng);
	void BasisFunction(double u, double v, int nen, const vector<array<double, phy_dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], double &detJ);
	void PointFormValue(vector<double> &Nx, const vector<array<double, dof_all>> &U, double Value[dof_all]);
	void PointFormGrad(vector<array<double, dim>> &dNdx, const vector<array<double, dof_all>> &U, double Value[dof_all][dim]);
	void PointFormHess(vector<array<array<double, dim>, dim>> &d2Ndx2, const vector<array<double, dof_all>> &U, double Value[dof_all][dim][dim]);
	void Tau(double J[dim][dim], double u[dof_all], double &tauM, double &tauC);
	void FineScale(double tauM, double tauC, double u[dim], double u_x[dim], double u_y[dim], double u_xx[dim], double u_yy[dim], double p, double p_x, double p_y, double u_s[dim], double &p_s);
	void Residual(vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], const double detJ, const vector<array<double, dof_all>> &U, vector<array<double, dof_all>> Re);
	void Tangent(vector<double> &Nx, vector<array<double, dim>> &dNdx, double dudx[dim][dim], const double detJ, const vector<array<double, dof_all>> &U, vector<array<vector<array<double, dof_all>>, dof_all>> &Ke);
	void BuildLinearSystemProcess(const vector<Vertex2D> &cpts, const vector<array<double, dim>> &velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, dof_all>>, dof_all>> &Ke, vector<array<double, dof_all>> &Re);
	void MatrixAssembly(vector<array<vector<array<double, dof_all>>, dof_all>> &Ke, const vector<int> &IEN, Mat &GK);
	void ResidualAssembly(vector<array<double, dof_all>> &Re, const vector<int> &IEN, Vec &GR);

	/*Postprocessing*/
	void ResultCal_Bezier(double u, double v, const Element2D &bzel, double pt[phy_dim], double result[dof_all], double dudx[dim], double &detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex2D> &pts, const vector<Element2D> &mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 4>> sele, int step, string fn);

public:
	/*Preprocessing*/
	void InitializeProblem(const int ndof, const int n_bz, const vector<double> &Vel0, const vector<double> &Pre0, const vector<double> &var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex2D> &cpts, const vector<Element2D> &tmesh, const vector<array<double, dim>> &velocity_bc, string fn);
};

#endif