#ifndef TRANSPORTOPT2D_H
#define TRANSPORTOPT2D_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "BasicDataStructure.h"
#include "time.h"

using namespace std;

const int degree = 3;
const int dim = 2;
const int state_num = 7;
const int ctrl_num = 4;

const int bzpt_num = 16;


class TransportOpt2D
{
private:
	vector<double> Gpt;
	vector<double> wght;
	const double PI = 4 * atan(1.0);

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
	SNES snes;
	PC pc;

	
	// Mat GK;
	// Vec GR;
	Vec temp_solution;
	
	double dt;
	double velocity_max;
	double nu;
	double rou;
	double alphaM;
	double alphaF;
	double Gama;
	double fx, fy, fz;//elemental body force


	// Scale of the problem
	int nTstep, nPoint;

	// constant parameters
	vector<double> par; //Dn0, v_plus, v_minus, k+, k-,k'+,k'-,c0,mu
	double alpha1, alpha2;
	double beta1, beta2;

	// state variables
	vector<double> n0, n_plus, n_minus;
	vector<double> Vel_plus[dim], Vel_minus[dim];
	// control variables
	vector<double> f_plus[dim], f_minus[dim];
	// penalty variables
	vector<double> lambda[3 + 2 * dim];

	// Unit Np * Np Matrix
	Mat M, K, P[dim];
	
	// KKT system matrix
	Mat Asubmat[9];

	// 
	// Vec Residual;
	// Mat Tangent;

	vector<double> Vel;
	vector<double> Pre;
	vector<double> gh; //prescribed boundary velocities

public:
	TransportOpt2D();
private:
	void BodyForce(double x, double y, double &Fx, double &Fy);
	void ReadBezierElementProcess(string fn);

	/*Analysis*/
	void GaussInfo(int ng);	
	void BasisFunction(double u, double v, int nen, const vector<array<double, dim>>& pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], double& detJ);

	// void BasisFunction(double u, double v, double w, int nen, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, vector<array<array<double, 3>, 3>> &dN2dx2, double dudx[3][3], double& detJ);
	void PointFormValue(vector<double> &Nx, const vector<double> &U, double Value);
	void PointFormGrad(vector<array<double, dim>>& dNdx, const vector<double> &U, double Value[dim]);
	void PointFormHess(vector<array<array<double, dim>, dim>>& d2Ndx2, const vector<double> &U, double Value[dim][dim]);

	void ComputeMassMatrix(vector<double>& Nx, const double detJ, vector<vector<double>>& MassMat);
	void ComputeStiffMatrix(vector<array<double, dim>>& dNdx, const double detJ, vector<vector<double>>& StiffMat);
	void ComputeParMatrix(vector<double>& Nx, vector<array<double, dim>>& dNdx, const double detJ, int dir, vector<vector<double>>& ParMat);

	void GetMatrixPosition(int row, int n_var, int &i_point, int &i_var, int &i_tstep);
	void FormMatrixA11(Mat M, Mat K, Mat &A);
	void FormMatrixA12(Mat M, Mat K, Mat &A);
	void FormMatrixA13(Mat M, Mat K, Mat P[dim], Mat &A);
	void FormMatrixA21(Mat M, Mat K, Mat &A);
	void FormMatrixA22(Mat M, Mat K, Mat &A);
	void FormMatrixA23(Mat M, Mat K, Mat &A);
	void FormMatrixA31(Mat M, Mat K, Mat P[dim], Mat &A);
	void FormMatrixA32(Mat M, Mat K, Mat &A);
	void FormMatrixA33(Mat M, Mat K, Mat &A);


	//void ComputeConvectionMatrix(vector<double>& Nx, vector<array<double, 3>>& dNdx, const double detJ, const vector<double> &U, vector<vector<double>> ConvectMat);
	void TangentAssemble();
	void ResidualAssemble();

	void Residual(vector<double>& Nx, vector<array<double, 3>>& dNdx, vector<array<array<double, 3>, 3>>& dN2dx2, double dudx[3][3], const double detJ, const vector<array<double, 4>> &U, vector<array<double, 4>> Re);
	void Tangent(vector<double> &Nx, vector<array<double, 3>>& dNdx, double dudx[3][3], const double detJ, const vector<array<double, 4>>& U, vector<array<vector<array<double, 4>>, 4>>& Ke);
	void BuildLinearSystemProcess(const vector<Vertex2D>& cpts, const vector<array<double, 2>>& velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, 4>>, 4>>& Ke, vector<array<double, 4>> &Re);

	void MatrixAssembly(vector<vector<double>>Emat, const vector<int>& IEN, Mat& Gmat);



	void TangentAssembly(vector<array<vector<array<double, 4>>, 4>>& Ke, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<array<double,4>> &Re, const vector<int>& IEN, Vec& GR);

	/*Postprocessing*/
	void ResultCal_Bezier(double u, double v, const Element2D& bzel, double pt[3], double result[4], double dudx[3], double& detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex2D>& pts, const vector<Element2D>& mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 8>> sele, int step, string fn);
public:
	/*Preprocessing*/
	void InitializeProblem(const int ndof, const int n_bz, const vector<double>& Vel0, const vector<double>& Pre0, const vector<double>& var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex2D>& cpts, const vector<Element2D>& tmesh, const vector<array<double, 2>>& velocity_bc, string fn);
};

#endif