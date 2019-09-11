#include "TransportOpt2D.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

double MatrixDet(double dxdt[2][2])
{
	double det = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
	return det;
}

double MatrixDet(double dxdt[3][3])
{
	double det = dxdt[0][0] * dxdt[1][1] * dxdt[2][2] + dxdt[0][1] * dxdt[1][2] * dxdt[2][0] + dxdt[0][2] * dxdt[2][1] * dxdt[1][0] -
		(dxdt[0][2] * dxdt[1][1] * dxdt[2][0] + dxdt[0][0] * dxdt[1][2] * dxdt[2][1] + dxdt[1][0] * dxdt[0][1] * dxdt[2][2]);
	
	return det;
}

void Matrix2DInverse(double dxdt[2][2], double dtdx[2][2])
{
	double det = MatrixDet(dxdt);
	dtdx[0][0] = 1 / det*(dxdt[1][1]);
	dtdx[0][1] = 1 / det*(-dxdt[0][1]);
	dtdx[1][0] = 1 / det*(-dxdt[1][0]);
	dtdx[1][1] = 1 / det*(dxdt[0][0]);
}

void Matrix3DInverse(double dxdt[3][3], double dtdx[3][3])
{
	double det = MatrixDet(dxdt);
	dtdx[0][0] = 1 / det*(dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]);
	dtdx[0][1] = 1 / det*(dxdt[2][1] * dxdt[0][2] - dxdt[0][1] * dxdt[2][2]);
	dtdx[0][2] = 1 / det*(dxdt[0][1] * dxdt[1][2] - dxdt[1][1] * dxdt[0][2]);
	dtdx[1][0] = 1 / det*(dxdt[2][0] * dxdt[1][2] - dxdt[1][0] * dxdt[2][2]);
	dtdx[1][1] = 1 / det*(dxdt[0][0] * dxdt[2][2] - dxdt[0][2] * dxdt[2][0]);
	dtdx[1][2] = 1 / det*(dxdt[1][0] * dxdt[0][2] - dxdt[0][0] * dxdt[1][2]);
	dtdx[2][0] = 1 / det*(dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	dtdx[2][1] = 1 / det*(dxdt[0][1] * dxdt[2][0] - dxdt[0][0] * dxdt[2][1]);
	dtdx[2][2] = 1 / det*(dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0]);
}

TransportOpt2D::TransportOpt2D()
{
	comm = MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;	
}

void TransportOpt2D::BodyForce(double x, double y, double &Fx, double &Fy)
{
	Fx = 0;
	Fy = 0;
}

void TransportOpt2D::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";

	cout << fname_cmat <<endl;

	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open())	{
		fin_cmat >> neles;
		cout << neles << endl;
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i<neles; i++){
			if (i == ele_process[add]){
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (int j = 0; j < nfunctions; j++){
					for (int k = 0; k < 16; k++){
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else{
				fin_cmat >> stmp >> nfunctions >> itmp;
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (int j = 0; j < nfunctions; j++)
					for (int k = 0; k < 16; k++)
						fin_cmat >> stmp;
			}
		}
		fin_cmat.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	cout << fname_bzpt <<endl;
	add = 0;
	if (fin_bzpt.is_open()){
		fin_bzpt >> npts;
		cout << npts << endl;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++)	{
			if (e == ele_process[add]){
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++){
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
				}
				add++;
			}
			else{
				for (int i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}
		}
		fin_bzpt.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

void TransportOpt2D::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
	case 2:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;							wght[1] = 1.;
		break;
	}
	case 3:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.1127016653792583;			Gpt[1] = 0.5;							Gpt[2] = 0.8872983346207417;
		wght[0] = 0.5555555555555556;			wght[1] = 0.8888888888888889;			wght[2] = 0.5555555555555556;
		break;
	}
	case 4:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.06943184420297371;			Gpt[1] = 0.33000947820757187;			Gpt[2] = 0.6699905217924281;			Gpt[3] = 0.9305681557970262;
		wght[0] = 0.3478548451374539;			wght[1] = 0.6521451548625461;			wght[2] = 0.6521451548625461;			wght[3] = 0.3478548451374539;
		break;
	}
	case 5:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;		Gpt[2] = 0.5;						Gpt[3] = 0.7692346550528415;		Gpt[4] = 0.953089922969332;
		wght[0] = 0.2369268850561891;		wght[1] = 0.4786286704993665;		wght[2] = 0.5688888888888889;		wght[3] = 0.4786286704993665;		wght[4] = 0.2369268850561891;
		break;
	}
	default:
	{
		Gpt.resize(2);
		wght.resize(2);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;							wght[1] = 1.;
		break;
	}
	}
}

void TransportOpt2D::InitializeProblem(const int ndof, const int n_bz, const vector<double>& Vel0, const vector<double>& Pre0, const vector<double>& var)
{
	
	MPI_Barrier(comm);

	PetscPrintf(PETSC_COMM_WORLD, "Initializing...\n");
	
	/*Initialize parameters*/
	GaussInfo(4);
	n_bzmesh = n_bz;

	// Scale of the problem
	nPoint = ndof;	

	// constant parameters
	alpha0 = var[9];
	alpha1 = var[10];
	alpha2 = var[11];
	beta1 = var[12];
	beta2 = var[13];
	dt = var[14];
	nTstep = var[15];
	par = var;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-

	// state variables
	n0.resize(nPoint * nTstep);
	n_plus.resize(nPoint * nTstep);
	n_minus.resize(nPoint * nTstep);
	for(int i = 0; i < dim; i++)
	{
		Vel_plus[i].resize(nPoint * nTstep);
		Vel_minus[i].resize(nPoint * nTstep);
	}
	// control variables
	for(int i = 0; i < dim; i++)
	{
		f_plus[i].resize(nPoint * nTstep);
		f_minus[i].resize(nPoint * nTstep);
	}

	// penalty variables
	for(int i = 0; i < 3 + 2 * dim; i++)
		lambda[i].resize(nPoint * nTstep);

	cout << "nPoint: "<< ndof << endl;
	cout << "nTstep: "<< nTstep << endl;

	/*Initialize petsc vector, matrix*/

	ierr = MatCreate(PETSC_COMM_WORLD, &M); 
	ierr = MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
	ierr = MatSetType(M, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(M, 100, NULL, 100, NULL);
	ierr = MatSetUp(M); 

	ierr = MatCreate(PETSC_COMM_WORLD, &K); 
	ierr = MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
	ierr = MatSetType(K, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(K, 100, NULL, 100, NULL);
	ierr = MatSetUp(K); 

	ierr = MatCreate(PETSC_COMM_WORLD, &P[0]); 
	ierr = MatSetSizes(P[0], PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
	ierr = MatSetType(P[0], MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(P[0], 100, NULL, 100, NULL);	
	ierr = MatSetUp(P[0]); 

	ierr = MatCreate(PETSC_COMM_WORLD, &P[1]); 
	ierr = MatSetSizes(P[1], PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
	ierr = MatSetType(P[1], MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(P[1], 100, NULL, 100, NULL);
	ierr = MatSetUp(P[1]); 

	// ierr = MatSetOption(GK, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	// ierr = MatSetUp(GK); 

	cout << "Setup Initial Vector\n";

	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &Y_k);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * ctrl_num * nTstep,  &U_k);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &L_k);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &Y_d);

	ierr = VecSet(Y_k, 1.0);
	ierr = VecSet(U_k, 1.0);
	ierr = VecSet(L_k, 1.0);
	ierr = VecSet(Y_d, 0.0);

	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &Res_nl);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * (2 * state_num + ctrl_num) * nTstep, &temp_solution);
	ierr = VecSet(Res_nl, 0.0);

	ierr = VecCreate(PETSC_COMM_WORLD, &ResVec);
	ierr = VecSetSizes(ResVec, PETSC_DECIDE, nPoint * nTstep * (state_num * 2 + ctrl_num));
	ierr = VecSetFromOptions(ResVec);


}

/// Need to pay attention that basis function right now is computed for iga using Bezier element definition
void TransportOpt2D::BasisFunction(double u, double v, int nen, const vector<array<double, dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], double &detJ)
{
	double Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
	double Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};

	double dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
	double dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};

	double dN2du2[4] = {6. * (1. - u), -12. + 18. * u, 6. - 18. * u, 6. * u};
	double dN2dv2[4] = {6. * (1. - v), -12. + 18. * v, 6. - 18. * v, 6. * v};

	double dNdt[bzpt_num][dim];
	double dN2dt2[bzpt_num][dim][dim];
	double Nx_bz[bzpt_num];
	double dNdx_bz[bzpt_num][dim];
	double dN2dx2_bz[bzpt_num][dim][dim];

	Nx.clear();
	dNdx.clear();
	dN2dx2.clear();
	Nx.resize(nen, 0);
	dNdx.resize(nen, {0});
	dN2dx2.resize(nen, {{0}});

	int i, j, k, a, b, c, loc;
	loc = 0;

	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < 4; k++)
		{
			Nx_bz[loc] = Nu[k] * Nv[j];
			dNdt[loc][0] = dNdu[k] * Nv[j];
			dNdt[loc][1] = Nu[k] * dNdv[j];
			dN2dt2[loc][0][0] = dN2du2[k] * Nv[j];
			dN2dt2[loc][0][1] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][0] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j];
			loc++;
		}
	}

	double dxdt[2][2] = {{0}};
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a < 2; a++)
			for (b = 0; b < 2; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

	double dtdx[2][2] = {{0}};
	Matrix2DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i < bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
	}
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.25 * detJ;

	//2nd derivatives
	double dx2dt2[2][4] = {{0}};
	double dt2dx2[2][4] = {{0}};
	for (int l = 0; l < 2; l++)
		for (loc = 0; loc < bzpt_num; loc++)
			for (a = 0; a < 2; a++)
				for (b = 0; b < 2; b++)
					dx2dt2[l][2 * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				for (a = 0; a < 2; a++)
					for (b = 0; b < 2; b++)
						for (c = 0; c < 2; c++)
							dt2dx2[c][2 * i + j] -= dx2dt2[k][2 * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

	for (loc = 0; loc < bzpt_num; loc++)
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				dN2dx2_bz[loc][i][j] = 0.;

	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 2; j++)
			{
				for (a = 0; a < 2; a++)
				{
					for (b = 0; b < 2; b++)
					{
						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
					}
					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][2 * i + j];
				}
			}
		}
	}

	for (i = 0; i < nen; i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 2; m++)
			{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
				for (int n = 0; n < 2; n++)
				{
					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
				}
			}
		}
	}
}

// For 3D
// void TransportOpt2D::BasisFunction(double u, double v, double w, int nen, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, vector<array<array<double, 3>, 3>> &dN2dx2, double dudx[3][3], double& detJ)
// {
// 	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
// 	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
// 	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };

// 	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
// 	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
// 	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	
// 	double dN2du2[4] = { 6.*(1. - u),-12. + 18.*u,6. - 18.*u,6.*u };
// 	double dN2dv2[4] = { 6.*(1. - v),-12. + 18.*v,6. - 18.*v,6.*v };
// 	double dN2dw2[4] = { 6.*(1. - w),-12. + 18.*w,6. - 18.*w,6.*w };

// 	double dNdt[bzpt_num][3];
// 	double dN2dt2[bzpt_num][3][3];
// 	double Nx_bz[bzpt_num];
// 	double dNdx_bz[bzpt_num][3];
// 	double dN2dx2_bz[bzpt_num][3][3];

// 	Nx.clear();
// 	dNdx.clear();
// 	dN2dx2.clear();
// 	Nx.resize(nen, 0);
// 	dNdx.resize(nen, { 0 });
// 	dN2dx2.resize(nen, { {0} });

// 	int i, j, k, a, b, c, loc;
// 	loc = 0;
// 	for (i = 0; i<4; i++){
// 		for (j = 0; j<4; j++){
// 			for (k = 0; k < 4; k++)	{
// 				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
// 				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
// 				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
// 				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
// 				dN2dt2[loc][0][0] = dN2du2[k] * Nv[j] * Nw[i];				dN2dt2[loc][0][1] = dNdu[k] * dNdv[j] * Nw[i];				dN2dt2[loc][0][2] = dNdu[k] * Nv[j] * dNdw[i];
// 				dN2dt2[loc][1][0] = dNdu[k] * dNdv[j] * Nw[i];				dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j] * Nw[i];				dN2dt2[loc][1][2] = Nu[k] * dNdv[j] * dNdw[i];
// 				dN2dt2[loc][2][0] = dNdu[k] * Nv[j] * dNdw[i];				dN2dt2[loc][2][1] = Nu[k] * dNdv[j] * dNdw[i];				dN2dt2[loc][2][2] = Nu[k] * Nv[j] * dN2dw2[i];
// 				loc++;
// 			}
// 		}
// 	}

// 	double dxdt[3][3] = { {0} };
// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (a = 0; a<3; a++)	
// 			for (b = 0; b<3; b++)
//  				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

// 	double dtdx[3][3] = { { 0 } };
// 	Matrix3DInverse(dxdt, dtdx);

// 	//1st derivatives
// 	for (i = 0; i<bzpt_num; i++){
// 		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0] + dNdt[i][2] * dtdx[2][0];
// 		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1] + dNdt[i][2] * dtdx[2][1];
// 		dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
// 	}
// 	for (int i = 0; i < 3; i++)
// 		for (int j = 0; j < 3; j++)
// 			dudx[i][j] = dtdx[i][j];

// 	detJ = MatrixDet(dxdt);
// 	detJ = 0.125*detJ;

// 	//2nd derivatives
// 	double dx2dt2[3][9] = { {0} };
// 	double dt2dx2[3][9] = { {0} };
// 	for (int l = 0; l < 3; l++)	
// 		for (loc = 0; loc<bzpt_num; loc++)
// 			for (a = 0; a<3; a++)
// 				for (b = 0; b<3; b++)
// 					dx2dt2[l][3*a+b] += pt[loc][l] * dN2dt2[loc][a][b];

// 	for (i = 0; i < 3; i++)
// 		for (j = 0; j < 3; j++)
// 			for (k = 0; k < 3; k++)
// 				for (a = 0; a < 3; a++)
// 					for (b = 0; b < 3; b++)
// 						for (c = 0; c < 3; c++)
// 							dt2dx2[c][3*i+j] -= dx2dt2[k][3*a+b]*dtdx[a][i]*dtdx[b][j]*dtdx[c][k];


	
// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (i = 0; i < 3; i++)
// 			for (j = 0; j < 3; j++)
// 				dN2dx2_bz[loc][i][j] = 0.;

// 	for (loc = 0; loc < bzpt_num; loc++){
// 		for (i = 0; i < 3; i++)	{
// 			for (j = 0; j < 3; j++)	{
// 				for (a = 0; a<3; a++){
// 					for (b = 0; b<3; b++){
// 						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i]*dtdx[b][j];
// 					}
// 					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][3*i+j];
// 				}
// 			}
// 		}
// 	}
	
// 	for (i = 0; i < nen; i++){
// 		for (j = 0; j < bzpt_num; j++)	{
// 			Nx[i] += cmat[i][j] * Nx_bz[j];
// 			for (int m = 0; m < 3; m++)	{
// 				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
// 				for (int n = 0; n < 3; n++)	{
// 					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
// 				}
// 			}
// 		}
// 	}
// }

void TransportOpt2D::PointFormValue(vector<double> &Nx, const vector<double> &U, double Value)
{
	Value = 0.;
	for (int j = 0; j < Nx.size(); j++)
		Value += U[j] * Nx[j];
}

void TransportOpt2D::PointFormGrad(vector<array<double, dim>>& dNdx, const vector<double> &U, double Value[dim])
{
	for (int j = 0; j < dim; j++)	
		Value[j] = 0.;

		for (int j = 0; j < dim; j++)	
			for (int k = 0; k < dNdx.size(); k++) 
				Value[j] += U[k] * dNdx[k][j];
}

void TransportOpt2D::PointFormHess(vector<array<array<double, dim>, dim>>& d2Ndx2, const vector<double> &U, double Value[dim][dim])
{

	for (int j = 0; j < dim; j++)	{
		for (int k = 0; k < dim; k++)	{
			Value[j][k] = 0.;
		}
	}

	for (int j = 0; j < dim; j++)	{
		for (int k = 0; k < dim; k++)	{
			for (int l = 0; l < d2Ndx2.size(); l++)	{
				Value[j][k] += U[l] * d2Ndx2[l][j][k];
			}
		}
	}

}

void TransportOpt2D::MatrixAssembly(vector<vector<double>>Emat, const vector<int>& IEN, Mat& Gmat)
{
	int i, j, A, B, m, n;	
	int row_start, row_end, row_now;
	int add=0;

	PetscInt *nodeList = new PetscInt[IEN.size()];
	PetscReal *tmpGmat = new PetscReal[IEN.size() * IEN.size()];
	
	for (m = 0; m<IEN.size(); m++){
		A = IEN[m];
		nodeList[m] = A;
		for (n = 0; n<IEN.size(); n++){
				B = IEN[n];
				tmpGmat[add] = Emat[m][n];
				add++;
		}
	}

	MatSetValues(Gmat, IEN.size(), nodeList, IEN.size(), nodeList, tmpGmat, ADD_VALUES);
	delete nodeList;
	delete tmpGmat;	
}

void TransportOpt2D::ComputeMassMatrix(vector<double>& Nx, const double detJ, vector<vector<double>>& MassMat)
{
	int a, b, nen = Nx.size();
	for (a = 0; a < nen; a++)
		for (b = 0; b < nen; b++)
			MassMat[a][b] = Nx[a] * Nx[b] * detJ;
}

void TransportOpt2D::ComputeStiffMatrix(vector<array<double, dim>>& dNdx, const double detJ, vector<vector<double>>& StiffMat)
{
	int a, b, c, nen = dNdx.size();
	double tmp = 0;
	for (a = 0; a < nen; a++)
	{
		for (b = 0; b < nen; b++)
		{
			tmp = 0;
			for (c = 0; c < dim; c++)
			{
				tmp += dNdx[a][c] * dNdx[b][c] * detJ;
			}
			StiffMat[a][b] = tmp;
		}
	}
}

void TransportOpt2D::ComputeParMatrix(vector<double>& Nx, vector<array<double, dim>>& dNdx, const double detJ, int dir, vector<vector<double>>& ParMat)
{
	int a, b, nen = Nx.size();
	for (a = 0; a < nen; a++)
		for (b = 0; b < nen; b++)
			ParMat[a][b] = Nx[a] * dNdx[b][dir] * detJ;
}

void TransportOpt2D::ComputeResVector(vector<double>& Nx, vector<array<double, dim>>& dNdx, const vector<int>& IEN, const double detJ)
{
	int a, b, nen = IEN.size();
	PetscInt *rows;
	PetscReal *vals;

	PetscMalloc1(IEN.size()*state_num, &rows);
	PetscMalloc1(IEN.size()*state_num, &vals);

	vector<double> n0_node, nplus_node, nminus_node;
	vector<double> vplus_node[dim], vminus_node[dim];
	n0_node.clear(); n0_node.resize(nen, 0);
	nplus_node.clear(); nplus_node.resize(nen, 0);
	nminus_node.clear(); nminus_node.resize(nen, 0);
	for(int i=0;i<dim;i++)
	{
		vplus_node[i].clear(); vplus_node[i].resize(nen, 0);
		vminus_node[i].clear(); vminus_node[i].resize(nen, 0);
	}

	for(int j = 0; j < nTstep; j++)
	{
		
		for(int i = 0; i < nen; i++)
		{
			int A = IEN[i] + j * nPoint;
			n0_node[i] = n0[A];
			nplus_node[i] = n_plus[A];
			nminus_node[i] = n_minus[A];
			vplus_node[0][i] = Vel_plus[0][A];
			vplus_node[1][i] = Vel_plus[1][A];
			vminus_node[0][i] = Vel_minus[0][A];
			vminus_node[1][i] = Vel_minus[1][A];
		}

		double n0Val;
		double npVal, npGrad[dim], vpxVal, vpxGrad[dim], vpyVal, vpyGrad[dim];
		double nmVal, nmGrad[dim], vmxVal, vmxGrad[dim], vmyVal, vmyGrad[dim];

		PointFormValue(Nx,n0_node,n0Val);

		PointFormValue(Nx,nplus_node,npVal);
		PointFormValue(Nx,vplus_node[0],vpxVal);
		PointFormValue(Nx,vplus_node[1],vpyVal);
		PointFormGrad(dNdx,nplus_node,npGrad);
		PointFormGrad(dNdx,vplus_node[0],vpxGrad);
		PointFormGrad(dNdx,vplus_node[1],vpyGrad);

		PointFormValue(Nx,nminus_node,nmVal);
		PointFormValue(Nx,vminus_node[0],vmxVal);
		PointFormValue(Nx,vminus_node[1],vmyVal);
		PointFormGrad(dNdx,nminus_node,nmGrad);
		PointFormGrad(dNdx,vminus_node[0],vmxGrad);
		PointFormGrad(dNdx,vminus_node[1],vmyGrad);

		for(int i = 0; i < nen; i++)
		{
			for(int k = 0; k < state_num; k++)
			{
				int A = IEN[i] + k * nPoint + j * nPoint * state_num;
				rows[i + nen * k] = A;
			}
			vals[i + nen * 0] = 0.;
			vals[i + nen * 1] = -Nx[i] * (vpxVal * npGrad[0] + vpyVal * npGrad[1]) * detJ * dt;
			vals[i + nen * 2] = -Nx[i] * (vmxVal * nmGrad[0] + vmyVal * nmGrad[1]) * detJ * dt;
			vals[i + nen * 3] = -Nx[i] * (vpxVal * vpxGrad[0] + vpyVal * vpxGrad[1]) * detJ * dt;
			vals[i + nen * 4] = -Nx[i] * (vpxVal * vpyGrad[0] + vpyVal * vpyGrad[1]) * detJ * dt;
			vals[i + nen * 5] = -Nx[i] * (vmxVal * vmxGrad[0] + vmyVal * vmxGrad[1]) * detJ * dt;
			vals[i + nen * 6] = -Nx[i] * (vmxVal * vmyGrad[0] + vmyVal * vmyGrad[1]) * detJ * dt;
		}
		
		if(j==0)
		{
			for(int i = 0; i < nen; i++)
			{
				vals[i + nen * 0] += Nx[i] * n0Val * detJ;
				vals[i + nen * 1] += Nx[i] * npVal * detJ;
				vals[i + nen * 2] += Nx[i] * nmVal * detJ;
				vals[i + nen * 3] += Nx[i] * vpxVal * detJ;
 				vals[i + nen * 4] += Nx[i] * vpyVal * detJ;
 				vals[i + nen * 5] += Nx[i] * vmxVal * detJ;
 				vals[i + nen * 6] += Nx[i] * vmyVal * detJ;
 			}
		}
		VecSetValues(Res_nl, nen * state_num, rows, vals, ADD_VALUES); 
	}

	PetscFree(rows);
	PetscFree(vals);

	delete rows, vals;
}


void TransportOpt2D::GetMatrixPosition(int row, int n_var, int &i_point, int &i_var, int &i_tstep)
{
    i_tstep = row / (n_var * nPoint);
    i_var = row % (n_var * nPoint) / nPoint;
    i_point = row % (n_var * nPoint) % nPoint;
}

void TransportOpt2D::FormMatrixA11(Mat M, Mat K, Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  
	PetscInt *col;
	PetscReal scale;
	PetscReal *vals;
	const PetscReal *vals_m, *vals_k;

	PetscMalloc1(nPoint, &col);
	PetscMalloc1(nPoint, &vals);
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	for(row = start; row < end; row++)
	{
        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
		MatGetRow(K,ind_point,NULL,NULL,&vals_k);
		MatGetRow(M,ind_point,NULL,NULL,&vals_m);
		
		switch (ind_var)
		{
		case 0:
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_k[i];
			scale = alpha0;
			break;
		case 1:
			// MatGetRow(K,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_k[i];
			scale = alpha1;
			break;
		case 2:
			// MatGetRow(K,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_k[i];
			scale = alpha2;
			break;
		case 3:
			// MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_m[i];
			scale = beta1;
			break;
		case 4:
			// MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_m[i];
			scale = beta1;
			break;
		case 5:
			// MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_m[i];
			scale = beta2;
			break;
		case 6:
			// MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
			for(int i = 0; i < nPoint; i++)
				vals[i] = vals_m[i];
			scale = beta2;
			break;
		// default:
		// 	for(int i = 0; i < nPoint; i++)
		// 		vals[i] = 0.;
		// 	break;
		}

		if(ind_time == 0 || ind_time == (nTstep - 1))
			scale = scale * 0.5 * dt;
		else
			scale = scale * dt;

		for(int i = 0; i < nPoint; i++)
		{
			vals[i] = vals[i]*scale;
			col[i] = i + ind_var * nPoint + ind_time * nPoint * state_num;
		}
			

		MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);

		// if(ind_var>=3)
			MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
		// else if(ind_var>0)
			MatRestoreRow(K,ind_point,NULL,NULL,&vals_k);

	}

	PetscFree(col);
	PetscFree(vals);

	delete col;
	delete vals;
	delete vals_m, vals_k;
	
}

void TransportOpt2D::FormMatrixA12(Mat M, Mat K, Mat &A)
{
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, 0, NULL, 0, NULL);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void TransportOpt2D::FormMatrixA13(Mat M, Mat K, Mat P[dim], Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  

	const PetscReal *vals_m, *vals_k, *vals_px, *vals_py;

	PetscInt *col;
	PetscReal *vals;
	PetscReal scale;
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	for(row = start; row < end; row++)
	{
        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
		MatGetRow(M,ind_point,NULL,NULL,&vals_m);
		MatGetRow(K,ind_point,NULL,NULL,&vals_k);
		MatGetRow(P[0],ind_point,NULL,NULL,&vals_px);
		MatGetRow(P[1],ind_point,NULL,NULL,&vals_py);
		if(ind_time > 0)
		{
			switch (ind_var)
			{
			case 0:
				// PetscInt *col = new PetscInt[4 * nPoint];
				// PetscReal *vals = new PetscReal[4 * nPoint];
				PetscMalloc1(4 * nPoint, &col);
				PetscMalloc1(4 * nPoint, &vals);
				for(int i = 0; i < nPoint; i++)
				{
					vals[i + 0 * nPoint] = -vals_m[i];
					vals[i + 1 * nPoint] = vals_m[i] * (1 + dt * par[3] + dt * par[4]);
					vals[i + 2 * nPoint] = vals_m[i] * (-dt * par[3]);
					vals[i + 3 * nPoint] = vals_m[i] * (-dt * par[4]);
					col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
					col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
					col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
					col[i + 3 * nPoint] = i + (ind_var + 2) * nPoint + ind_time * nPoint * state_num; 
				}
				MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
				PetscFree(col);
				PetscFree(vals);
				break;
			case 1:
				// PetscInt *col = new PetscInt[5 * nPoint];
				// PetscReal *vals = new PetscReal[5 * nPoint];
				PetscMalloc1(5 * nPoint, &col);
				PetscMalloc1(5 * nPoint, &vals);
				for(int i = 0; i < nPoint; i++)
				{
					vals[i + 0 * nPoint] = -vals_m[i];
					vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[5]);
					vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[5]);
					vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
					vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
					col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
					col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
					col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
					col[i + 3 * nPoint] = i + (ind_var + 3) * nPoint + ind_time * nPoint * state_num;
					col[i + 4 * nPoint] = i + (ind_var + 4) * nPoint + ind_time * nPoint * state_num; 
				}
				MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
				PetscFree(col);
				PetscFree(vals);
				break;
			case 2:
				// PetscInt *col = new PetscInt[5 * nPoint];
				// PetscReal *vals = new PetscReal[5 * nPoint];
				PetscMalloc1(5 * nPoint, &col);
				PetscMalloc1(5 * nPoint, &vals);
				for(int i = 0; i < nPoint; i++)
				{
					vals[i + 0 * nPoint] = -vals_m[i];
					vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[6]);
					vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[6]);
					vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
					vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
					col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
					col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
					col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
					col[i + 3 * nPoint] = i + (ind_var + 5) * nPoint + ind_time * nPoint * state_num;
					col[i + 4 * nPoint] = i + (ind_var + 6) * nPoint + ind_time * nPoint * state_num; 
				}
				MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
				PetscFree(col);
				PetscFree(vals);
				break;
			default:
				// PetscInt *col = new PetscInt[2 * nPoint];
				// PetscReal *vals = new PetscReal[2 * nPoint];
				PetscMalloc1(2 * nPoint, &col);
				PetscMalloc1(2 * nPoint, &vals);
				for(int i = 0; i < nPoint; i++)
				{
					vals[i + 0 * nPoint] = -vals_m[i];
					vals[i + 1 * nPoint] = vals_m[i] + dt * par[8] * vals_k[i];
					col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
					col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
				}
				MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
				PetscFree(col);
				PetscFree(vals);
				break;
			}
		}
		MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
		MatRestoreRow(K,ind_point,NULL,NULL,&vals_k);
		MatRestoreRow(P[0],ind_point,NULL,NULL,&vals_px);
		MatRestoreRow(P[1],ind_point,NULL,NULL,&vals_py);
	}
	delete col, vals;
	delete vals_m, vals_k, vals_px, vals_py;
}

void TransportOpt2D::FormMatrixA21(Mat M, Mat K, Mat &A)
{
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, ctrl_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, 0, NULL, 0, NULL);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void TransportOpt2D::FormMatrixA22(Mat M, Mat K, Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  
	PetscInt *col;

	PetscReal scale;
	PetscReal *vals;
	const PetscReal *vals_m;

	PetscMalloc1(nPoint, &col);
	PetscMalloc1(nPoint, &vals);
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, ctrl_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	for(row = start; row < end; row++)
	{
        GetMatrixPosition(row, ctrl_num, ind_point, ind_var, ind_time);
		MatGetRow(M,ind_point,NULL,NULL,&vals_m);
		switch (ind_var)
		{
		case 0:
			scale = beta1;
			break;
		case 1:
			scale = beta1;
			break;
		case 2:
			scale = beta2;
			break;
		case 3:
			scale = beta2;
			break;
		default:
			break;
		}

		if(ind_time == 0 || ind_time == nTstep - 1 )
			scale = scale * 0.5 * dt;
		else
			scale = scale * dt;

		for(int i = 0; i < nPoint; i++)
		{
			vals[i] = vals_m[i]*scale;
			col[i] = i + ind_var * nPoint + ind_time * nPoint * ctrl_num;
		}
			
		MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
		MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
	}

	PetscFree(col);
	PetscFree(vals);

	delete col;
	delete vals;
	delete vals_m;
}

void TransportOpt2D::FormMatrixA23(Mat M, Mat K, Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  
	PetscInt *col;

	PetscReal scale;
	PetscReal *vals;
	const PetscReal *vals_m;

	PetscMalloc1(nPoint, &col);
	PetscMalloc1(nPoint, &vals);
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, ctrl_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	for(row = start; row < end; row++)
	{
        GetMatrixPosition(row, ctrl_num, ind_point, ind_var, ind_time);
		MatGetRow(M,ind_point,NULL,NULL,&vals_m);

		scale = -dt;
		for(int i = 0; i < nPoint; i++)
		{
			vals[i] = vals_m[i]*scale;
			col[i] = i + (ind_var + 3) * nPoint + ind_time * nPoint * state_num;
		}
			
		MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);

		MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
	}

	PetscFree(col);
	PetscFree(vals);

	delete col;
	delete vals;
	delete vals_m;
}

void TransportOpt2D::FormMatrixA31(Mat M, Mat K, Mat P[dim], Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  
	PetscReal *vals_m = new PetscReal[nPoint];
	PetscReal *vals_k = new PetscReal[nPoint];
	PetscReal *vals_px = new PetscReal[nPoint];
	PetscReal *vals_py = new PetscReal[nPoint];
	PetscReal scale;
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	// for(row = start; row < end; row++)
	// {
    //     GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
	// 	MatGetRow(M,ind_point,NULL,NULL,&vals_m);
	// 	MatGetRow(K,ind_point,NULL,NULL,&vals_k);
	// 	MatGetRow(P[0],ind_point,NULL,NULL,&vals_px);
	// 	MatGetRow(P[1],ind_point,NULL,NULL,&vals_py);
	// 	if(ind_time > 0)
	// 	{
	// 		switch (ind_var)
	// 		{
	// 		case 0:
	// 			PetscInt *col = new PetscInt[4 * nPoint];
	// 			PetscReal *vals = new PetscReal[4 * nPoint];
	// 			for(int i = 0; i < nPoint; i++)
	// 			{
	// 				vals[i + 0 * nPoint] = -vals_m[i];
	// 				vals[i + 1 * nPoint] = vals_m[i] * (1 + dt * par[3] + dt * par[4]);
	// 				vals[i + 2 * nPoint] = vals_m[i] * (-dt * par[3]);
	// 				vals[i + 3 * nPoint] = vals_m[i] * (-dt * par[4]);
	// 				col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
	// 				col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 3 * nPoint] = i + (ind_var + 2) * nPoint + ind_time * nPoint * state_num; 
	// 			}
	// 			MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
	// 			delete col, vals;
	// 			break;
	// 		case 1:
	// 			PetscInt *col = new PetscInt[5 * nPoint];
	// 			PetscReal *vals = new PetscReal[5 * nPoint];
	// 			for(int i = 0; i < nPoint; i++)
	// 			{
	// 				vals[i + 0 * nPoint] = -vals_m[i];
	// 				vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[5]);
	// 				vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[5]);
	// 				vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
	// 				vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
	// 				col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
	// 				col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 3 * nPoint] = i + (ind_var + 3) * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 4 * nPoint] = i + (ind_var + 4) * nPoint + ind_time * nPoint * state_num; 
	// 			}
	// 			MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
	// 			delete col, vals;
	// 			break;
	// 		case 2:
	// 			PetscInt *col = new PetscInt[5 * nPoint];
	// 			PetscReal *vals = new PetscReal[5 * nPoint];
	// 			for(int i = 0; i < nPoint; i++)
	// 			{
	// 				vals[i + 0 * nPoint] = -vals_m[i];
	// 				vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[6]);
	// 				vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[6]);
	// 				vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
	// 				vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
	// 				col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
	// 				col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 2 * nPoint] = i + (ind_var + 1) * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 3 * nPoint] = i + (ind_var + 5) * nPoint + ind_time * nPoint * state_num;
	// 				col[i + 4 * nPoint] = i + (ind_var + 6) * nPoint + ind_time * nPoint * state_num; 
	// 			}
	// 			MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
	// 			delete col, vals;
	// 			break;
	// 		default:
	// 			PetscInt *col = new PetscInt[2 * nPoint];
	// 			PetscReal *vals = new PetscReal[2 * nPoint];
	// 			for(int i = 0; i < nPoint; i++)
	// 			{
	// 				vals[i + 0 * nPoint] = -vals_m[i];
	// 				vals[i + 1 * nPoint] = vals_m[i] + dt * par[8] * vals_k[i];
	// 				col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) * nPoint * state_num;
	// 				col[i + 1 * nPoint] = i + ind_var * nPoint + ind_time * nPoint * state_num;
	// 			}
	// 			MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
	// 			delete col, vals;
	// 			break;
	// 		}
	// 	}
	// }

	delete vals_m, vals_k, vals_px, vals_py;
}

void TransportOpt2D::FormMatrixA32(Mat M, Mat K, Mat &A)
{
	PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;  
	PetscInt *col = new PetscInt[nPoint];
	PetscReal *vals = new PetscReal[nPoint];
	const PetscReal *vals_m;
	PetscReal scale;
	
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, nPoint, NULL, nPoint, NULL);
	MatGetOwnershipRange(A,&start,&end);
	
	for(row = start; row < end; row++)
	{
        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
		if(ind_var >= 3)
		{
			MatGetRow(M,ind_point,NULL,NULL,&vals_m);
			scale = -dt;
			for(int i = 0; i < nPoint; i++)
			{
				vals[i] = vals_m[i]*scale;
				col[i] = i + (ind_var - 3) * nPoint + ind_time * nPoint * ctrl_num;
			}
			MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
			MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
		}
	}

	delete col;
	delete vals;
}

void TransportOpt2D::FormMatrixA33(Mat M, Mat K, Mat &A)
{
	MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, 0, NULL, 0, NULL);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void TransportOpt2D::TangentMatSetup()
{
	PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;

	/// A11
	FormMatrixA11(M,K,Asubmat[0]);  // need debug
	MatAssemblyBegin(Asubmat[0], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[0], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA11.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[0],viewer);

	// A12 & A21 // No problem
	FormMatrixA12(M,K,Asubmat[1]);  // need debug
	MatAssemblyBegin(Asubmat[1], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[1], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA12.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[1],viewer);
	MatTranspose(Asubmat[1],MAT_INITIAL_MATRIX,&Asubmat[3]); //correct
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA21.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[3],viewer);

	// A13 & A31 // No problem
	FormMatrixA13(M,K,P,Asubmat[2]);  // need debug
	MatAssemblyBegin(Asubmat[2], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[2], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA13.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[2],viewer);
	MatTranspose(Asubmat[2],MAT_INITIAL_MATRIX,&Asubmat[6]); //correct
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA31.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[6],viewer);

	// A22  // No problem
	FormMatrixA22(M,K,Asubmat[4]);  
	MatAssemblyBegin(Asubmat[4], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[4], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA22.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[4],viewer);

	// A23 & A32 // No problem
	FormMatrixA23(M,K,Asubmat[5]); 
	MatAssemblyBegin(Asubmat[5], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[5], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA23.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[5],viewer);
	MatTranspose(Asubmat[5],MAT_INITIAL_MATRIX,&Asubmat[7]);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA32.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[7],viewer);

	// A33 // No problem
	FormMatrixA33(M,K,Asubmat[8]);  // need debug
	MatAssemblyBegin(Asubmat[8], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Asubmat[8], MAT_FINAL_ASSEMBLY);
	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matA33.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(Asubmat[8],viewer);



	MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, Asubmat, &TanMat);

	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/TanMat.output",&viewer);
    // ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(TanMat,viewer);

    PetscPrintf(PETSC_COMM_WORLD,  "finish Tangent Mat\n");

}

void TransportOpt2D::FormResVecb1()
{
	Vec vec_tmp0, vec_tmp1, vec_tmp2;
	VecDuplicate(Y_d, &bsub[0]);
	VecDuplicate(Y_d, &vec_tmp0);
	VecDuplicate(Y_d, &vec_tmp1);
	VecDuplicate(Y_d, &vec_tmp2);
	// VecCreateMPI();

	VecWAXPY(vec_tmp0, -1.0, Y_k, Y_d); //y_d - y_k
	MatMult(Asubmat[0], vec_tmp0, vec_tmp1); //vec_tmp1 = A11 * (y_d-y_k)

	MatMult(Asubmat[2], L_k, vec_tmp2); // vec_tmp2 = A13 * l_k
	VecWAXPY(bsub[0], -1.0, vec_tmp2, vec_tmp1); //vec_tmp1 - vec_tmp2 = A11 * (y_d-y_k) - A13 * l_k

	PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/VecB1.output",&viewer);
    ierr = PetscViewerPushFormat(viewer,format);
    ierr = VecView(bsub[0],viewer);

}
void TransportOpt2D::FormResVecb2()
{
	Vec vec_tmp0, vec_tmp1, vec_tmp2;

	VecDuplicate(U_k, &bsub[1]);
	VecDuplicate(U_k, &vec_tmp0);
	VecDuplicate(U_k, &vec_tmp1);
	VecDuplicate(U_k, &vec_tmp2);

	MatMult(Asubmat[4], U_k, vec_tmp1); //A22 * u_k
	VecScale(vec_tmp1, -1.0); // vec_tmp1 = -A22 * u_k

	MatMult(Asubmat[5], L_k, vec_tmp2); // vec_tmp2 = A13 * l_k
	VecWAXPY(bsub[1], -1.0, vec_tmp1, vec_tmp0); //vec_tmp1 - vec_tmp2 = -A22 * u_k - A23 * l_k

	PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/VecB2.output",&viewer);
    ierr = PetscViewerPushFormat(viewer,format);
    ierr = VecView(bsub[1],viewer);
}
void TransportOpt2D::FormResVecb3()
{
	Vec vec_tmp0, vec_tmp1, vec_tmp2, vec_tmp3;
	VecDuplicate(L_k, &vec_tmp0);
	VecDuplicate(L_k, &vec_tmp1);
	VecDuplicate(L_k, &vec_tmp2);
	VecDuplicate(L_k, &vec_tmp3);
	VecDuplicate(L_k, &bsub[2]);

	MatMult(Asubmat[6], Y_k, vec_tmp1); //A31 * y_k
	VecScale(vec_tmp1, -1.0); // vec_tmp1 = -A31 * y_k

	MatMult(Asubmat[7], U_k, vec_tmp2); //A32 * u_k
	VecScale(vec_tmp2, -1.0); // vec_tmp2 = -A32 * u_k

	VecCopy(Res_nl, vec_tmp3);
	VecWAXPY(vec_tmp0, 1.0, vec_tmp1, vec_tmp2); //vec_tmp1 + vec_tmp2 = -A31 * y_k -A32 * u_k
	VecWAXPY(bsub[2], 1.0, vec_tmp0, vec_tmp3); //vec_tmp0 + vec_tmp3 = d(y_k) -A31 * y_k -A32 * u_k

	PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/VecB3.output",&viewer);
    ierr = PetscViewerPushFormat(viewer,format);
    ierr = VecView(bsub[2],viewer);

}

void TransportOpt2D::ResidualVecSetup()
{
    FormResVecb1();
	FormResVecb2();
	FormResVecb3();

	PetscInt low, high, ldim, iglobal, i, k;
	PetscReal val, *array;
	PetscInt shift[3] = {0, state_num * nTstep * nPoint, (state_num + ctrl_num) * nTstep * nPoint};


	// VecGetOwnershipRange(bsub[0], &low, &high);
	// VecGetLocalSize(bsub[0], &ldim);
	// VecGetArray(bsub[0], &array);

	// for(i = 0; i < ldim; i ++)
	// {
	// 	val  = array[i];
	// 	iglobal = i + low;
	// 	VecSetValues(ResVec, 1, &iglobal, &val, INSERT_VALUES);
	// }

	for (k = 0; k < 3; k++)
	{
		VecGetOwnershipRange(bsub[k], &low, &high);
		VecGetLocalSize(bsub[k], &ldim);
		VecGetArray(bsub[k], &array);

		for (i = 0; i < ldim; i++)
		{
			val = array[i];
			iglobal = i + low + shift[k];
			VecSetValues(ResVec, 1, &iglobal, &val, INSERT_VALUES);
		}
	}
	VecAssemblyBegin(ResVec);
	VecAssemblyEnd(ResVec);

	// VecCreate(PETSC_COMM_WORLD, &ResVec);
	// VecSetSizes(ResVec, PETSC_DECIDE, nPoint * nTstep * (state_num * 2 + ctrl_num));
	// VecSetFromOptions(ResVec);
	
	// VecCreateNest(PETSC_COMM_WORLD, 3, NULL, bsub, &ResVec);
	// PetscPrintf(PETSC_COMM_WORLD,  "finish Residual Vec\n");

	PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/ResVec.output",&viewer);
    ierr = PetscViewerPushFormat(viewer,format);
    ierr = VecView(ResVec,viewer);

    PetscPrintf(PETSC_COMM_WORLD,  "finish Tangent Mat\n");

}

void TransportOpt2D::BuildLinearSystemProcess(const vector<Vertex2D>& cpts, const vector<array<double, 2>>& velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node)
{
	/*Build linear system in each process*/
	int e;
	cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "<< bzmesh_process.size() <<" elements.\n";
	for (e=0;e<bzmesh_process.size();e++){
		int nen, A;
		double ux_bc, uy_bc, uz_bc, p_bc;
	
		double dudx[dim][dim];
		double detJ;
		vector<double> Nx;
		vector<array<double, 2>> dNdx;
		vector<array<array<double, 2>, 2>> dN2dx2;
		vector<vector<double>> Mtmp, Ktmp, Pxtmp, Pytmp;

		// vector<array<double, 4>> Re;
		// vector<array<vector<array<double, 4>>, 4>> Ke;
	
		// vector<array<double, 4>> v_node;	

		vector<double> n0_node, nplus_node, nminus_node;
		vector<double> vplus_node[dim], vminus_node[dim];
		vector<double> fplus_node[dim], fminus_node[dim];
		
		nen = bzmesh_process[e].IEN.size();
		
		Nx.clear(); Nx.resize(nen, 0);
		dNdx.clear(); dNdx.resize(nen, { 0 });
		dN2dx2.clear(); dN2dx2.resize(nen, { {0} });

		n0_node.clear(); n0_node.resize(nen, 0);
		nplus_node.clear(); nplus_node.resize(nen, 0);
		nminus_node.clear(); nminus_node.resize(nen, 0);
		for(int i=0;i<dim;i++)
		{
			vplus_node[i].clear(); vplus_node[i].resize(nen, 0);
			vminus_node[i].clear(); vminus_node[i].resize(nen, 0);
			fplus_node[i].clear(); fplus_node[i].resize(nen, 0);
			fminus_node[i].clear(); fminus_node[i].resize(nen, 0);
		}
		
		Mtmp.clear(); Mtmp.resize(nen);
		Ktmp.clear(); Ktmp.resize(nen);
		Pxtmp.clear(); Pxtmp.resize(nen);
		Pytmp.clear(); Pytmp.resize(nen);

		for(int i = 0; i < nen; i++)
		{
			Mtmp[i].resize(nen, 0.);
			Ktmp[i].resize(nen, 0.);
			Pxtmp[i].resize(nen, 0.);
			Pytmp[i].resize(nen, 0.);
		}
	
		for (int i = 0; i < Gpt.size(); i++){
			for (int j = 0; j < Gpt.size(); j++){
					BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
					detJ = wght[i] * wght[j] * detJ;
					ComputeMassMatrix(Nx, detJ, Mtmp);
					ComputeStiffMatrix(dNdx, detJ, Ktmp);
					ComputeParMatrix(Nx, dNdx, detJ, 0, Pxtmp);
					ComputeParMatrix(Nx, dNdx, detJ, 1, Pytmp);

					ComputeResVector(Nx, dNdx, bzmesh_process[e].IEN, detJ);
					// Tangent(Nx, dNdx, dudx, detJ, v_node, Ke);
					// Residual(Nx, dNdx, dN2dx2, dudx, detJ, v_node, Re);
			}
		}

		/*Start element vector assembly*/
		
	
		/*Start element matrix assembly*/
		MatrixAssembly(Mtmp, bzmesh_process[e].IEN, M);
		MatrixAssembly(Ktmp, bzmesh_process[e].IEN, K);
		MatrixAssembly(Pxtmp, bzmesh_process[e].IEN, P[0]);
		MatrixAssembly(Pytmp, bzmesh_process[e].IEN, P[1]);
		// TangentAssembly(Ke, bzmesh_process[e].IEN, GK);
	
		// /*Apply Boundary Condition*/
		// for (int i = 0; i < nen; i++){
		// 	A = bzmesh_process[e].IEN[i];
		// 	if (cpts[A].label == 1) {
		// 		//inlet
		// 		ux_bc = velocity_max*velocity_bc[A][0] - velocity_node[A * 3];				
		// 		ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
		// 		uy_bc = velocity_max*velocity_bc[A][1] - velocity_node[A * 3 + 1];
		// 		ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
		// 		uz_bc = velocity_max*velocity_bc[A][2] - velocity_node[A * 3 + 2];
		// 		ApplyBoundaryCondition(uz_bc, i, 2, Ke, Re);
		// 	}
		// 	if (cpts[A].label == 0) {
		// 		//wall
		// 		ux_bc = 0.0 - velocity_node[A * 3];
		// 		ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
		// 		uy_bc = 0.0 - velocity_node[A * 3 + 1];
		// 		ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
		// 		uz_bc = 0.0 - velocity_node[A * 3 + 2];
		// 		ApplyBoundaryCondition(uz_bc, i, 2, Ke, Re);
		// 	}	
		// }

		// /*Start element vector assembly*/
		// ResidualAssembly(Re, bzmesh_process[e].IEN, GR);
	
		// /*Start element matrix assembly*/
		// TangentAssembly(Ke, bzmesh_process[e].IEN, GK);
	}
	cout << "Process " << comRank << " :complete build matrix and vector!\n";

	MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(P[0], MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(P[1], MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(Res_nl);
	// VecAssemblyBegin(GR);
	// MatAssemblyBegin(GK, MAT_FINAL_ASSEMBLY);	
}

void TransportOpt2D::ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, 4>>, 4>>& Ke, vector<array<double, 4>> &Re)
{
	int j, k;
	for (j = 0; j < Re.size(); j++)	{
		for (k = 0; k < 4; k++)		{
			Re[j][k] += bc_value*Ke[j][k][pt_num][variable_num];
		}
	}
	for (j = 0; j < Re.size(); j++)	{
		for (k = 0; k < 4; k++)		{
			Ke[pt_num][variable_num][j][k] = 0.0;		Ke[j][k][pt_num][variable_num] = 0.0;
		}
	}
	Re[pt_num][variable_num] = -bc_value;
	Ke[pt_num][variable_num][pt_num][variable_num] = 1.0;
}

void TransportOpt2D::AssignProcessor(vector<vector<int>> &ele_proc)
{
	/*Assign the partitioned bezier elements to this processor */
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void TransportOpt2D::Run(const vector<Vertex2D>& cpts, const vector<Element2D>& tmesh, const vector<array<double, 2>>& velocity_bc, string fn)
{
	
	int n_iterator(1);
	int l(0);
	time_t t0, t1;	
	vector<double> V_delta(3 * cpts.size()), P_delta(cpts.size());
	
	ReadBezierElementProcess(fn);
	for (l = 0; l<n_iterator; l++){
	
		PetscPrintf(PETSC_COMM_WORLD, "Iteration Step: %d\n", l);
		/*Build Linear System*/
		PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");	
		MatZeroEntries(M);
		MatZeroEntries(K);
		MatZeroEntries(P[0]);
		MatZeroEntries(P[1]);
		VecSet(Res_nl, 0.0);
		
		// VecSet(GR, 0);		
		t0 = time(NULL);
		BuildLinearSystemProcess(cpts, velocity_bc, Vel, Pre);
		// // VecAssemblyEnd(GR);
		// PetscPrintf(PETSC_COMM_WORLD, "Done Vector Assembly...\n");
		MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(P[0], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(P[1], MAT_FINAL_ASSEMBLY);
		VecAssemblyEnd(Res_nl);

		// PetscViewer viewer;
    	// PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
    	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matM.output",&viewer);
    	// ierr = PetscViewerPushFormat(viewer,format);
    	// ierr = MatView(M,viewer);
		// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matK.output",&viewer);
    	// ierr = PetscViewerPushFormat(viewer,format);
    	// ierr = MatView(K,viewer);
		// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matPx.output",&viewer);
    	// ierr = PetscViewerPushFormat(viewer,format);
    	// ierr = MatView(P[0],viewer);
		// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/matPy.output",&viewer);
    	// ierr = PetscViewerPushFormat(viewer,format);
    	// ierr = MatView(P[1],viewer);
		PetscPrintf(PETSC_COMM_WORLD, "Done Unit Matrix Assembly...\n");

		TangentMatSetup();
		ResidualVecSetup();
		
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 - t0);		
		MPI_Barrier(comm);
		PetscPrintf(PETSC_COMM_WORLD, "Solving...\n");

		t0 = time(NULL);
		
		/*Petsc KSP solver setting*/
		KSPCreate(PETSC_COMM_WORLD, &ksp);
		KSPSetOperators(ksp, TanMat, TanMat);
		// KSPGetPC(ksp, &pc);
		// {
		// 	MatNestGetISs(TanMat, isg, NULL);
		// 	PCFieldSplitSetIS(pc, "y", isg[0]);
		// 	PCFieldSplitSetIS(pc, "u", isg[1]);
		// 	PCFieldSplitSetIS(pc, "l", isg[2]);		

		// 	// PetscViewer viewer;
    	// 	// PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
    	// 	// ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/isg1.output",&viewer);
    	// 	// ierr = PetscViewerPushFormat(viewer,format);
		// 	// //  ISView(isg[0],viewer);
		// 	// // ISView(isg[0],iewer);
		// 	//   ISView(isg[2],viewer);
		// 	// //    ISView(isg[2],PETSC_VIEWER_STDOUT_WORLD);
		// }
		// PCSetType(pc, PCFIELDSPLIT);
		// PCSetUp(pc);
				
		
		// KSPSetType(ksp, KSPGMRES);

		// KSP *subksp;
		// PC subpc;
		// PetscInt first,nlocal;
		// KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, 100000);
		// // KSPSetPCSide(ksp, PC_RIGHT);
		// // KSPSetFromOptions(ksp);
		// KSPSetType(ksp, KSPGMRES);
		// KSPGMRESSetRestart(ksp, 500);
		// KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
		KSPSetUp(ksp);
		// PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
		// for (int i = 0; i<nlocal; i++) {
		// 	KSPGetPC(subksp[i], &subpc);		
		// 	PCSetType(subpc, PCILU);
		// 	KSPSetType(subksp[i], KSPGMRES);
		// 	KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE);
		// 	KSPSetPCSide(subksp[i], PC_RIGHT);
		// }
	
		/*Solving the equation*/

		// Vec ResVec_test;
		// ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * (2 * state_num + ctrl_num) * nTstep, &ResVec_test);
		// ierr = VecSet(ResVec_test, 1.0);

		// KSPSolve(ksp, ResVec_test, temp_solution);

		KSPSolve(ksp, ResVec, temp_solution);

		KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "iterations %d\n", its);
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
		PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Solving with time: %d \n", t1 - t0);

		PetscViewer viewer;
    	PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./debug/temp_solution.output",&viewer);
    	ierr = PetscViewerPushFormat(viewer,format);
    	ierr = VecView(temp_solution,viewer);

		// /*Collect the solution from all processors*/
		// Vec temp_solution_seq;
		// VecScatter ctx;
		// PetscReal    *_a;
		// VecScatterCreateToAll(temp_solution, &ctx, &temp_solution_seq);
		// VecScatterBegin(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		// VecScatterEnd(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		// VecGetArray(temp_solution_seq, &_a);		
		// MPI_Barrier(comm);
		
		// for (uint i = 0; i < cpts.size(); i++){
		// 	V_delta[3 * i] = PetscRealPart(_a[4 * i]);
		// 	V_delta[3 * i + 1] = PetscRealPart(_a[4 * i + 1]);
		// 	V_delta[3 * i + 2] = PetscRealPart(_a[4 * i + 2]);
		// 	P_delta[i] = PetscRealPart(_a[4 * i + 3]);
		// }
		// VecRestoreArray(temp_solution_seq, &_a);
		// VecScatterDestroy(&ctx);
		// VecDestroy(&temp_solution_seq);
		
		// for (uint i = 0; i < cpts.size(); i++){
		// 	Vel[3 * i] += V_delta[3 * i];
		// 	Vel[3 * i + 1] += V_delta[3 * i + 1];
		// 	Vel[3 * i + 2] += V_delta[3 * i + 2];
		// 	Pre[i] += P_delta[i];
		// }
		// MPI_Barrier(comm);
		// /*Visualize the result*/
		// if (comRank == 0){
		// 	cout << "Visualizing...\n";
		// 	VisualizeVTK_ControlMesh(cpts, tmesh, l, fn);
		// }
		// MPI_Barrier(comm);
		// VisualizeVTK_PhysicalDomain(l, fn + "final_physics");
	}
	// MatDestroy(&GK);
	// VecDestroy(&GR);
	VecDestroy(&temp_solution);
	KSPDestroy(&ksp);

	MPI_Barrier(comm);
}

void TransportOpt2D::VisualizeVTK_ControlMesh(const vector<Vertex2D>& spt, const vector<Element2D>& mesh, int step, string fn)
{
	stringstream ss;
	ss << step;
	
	//string fname = fn + "controlmesh_VelocityPressure_"+ss.str()+".vtk";
	string fname = fn + "controlmesh_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
				<< " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << Vel.size() / 3 << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i<Vel.size() / 3; i++)
		{
			fout << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 + 2] << "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<Pre.size(); i++)
		{
			fout << sqrt(Vel[i * 3] * Vel[i * 3] + Vel[i * 3 + 1] * Vel[i * 3 + 1] + Vel[i * 3 + 2] * Vel[i * 3 + 2]) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<Pre.size(); i++)
		{
			fout << Pre[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1 = fn + "velocityfield.txt";
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		for (uint i = 0; i<Vel.size() / 3; i++)
		{
			fout1 << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 + 2] << "\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void TransportOpt2D::VisualizeVTK_PhysicalDomain(int step, string fn)
{
	vector<array<double, 3>> spt_all;//sample points
	vector<double> sresult_all;
	vector<array<int, 8>> sele_all;
	double detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	double spt_proc[num_bzmesh_ele * 24];
	double sresult_proc[num_bzmesh_ele * 32];
	int sele_proc[num_bzmesh_ele * 8];
	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{

					double pt1[3], dudx[3];
					double result[4];
					ResultCal_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
					spt_proc[24 * e + loc*3 + 0] = pt1[0];
					spt_proc[24 * e + loc*3 + 1] = pt1[1];
					spt_proc[24 * e + loc*3 + 2] = pt1[2];
					sresult_proc[32 * e + loc * 4 + 0] = result[0];
					sresult_proc[32 * e + loc * 4 + 1] = result[1];
					sresult_proc[32 * e + loc * 4 + 2] = result[2];
					sresult_proc[32 * e + loc * 4 + 3] = result[3];
					loc++;
				
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					sele_proc[8 * e + 0] = 8 * e + a*nns[1] + b*ns + c;
					sele_proc[8 * e + 1] = 8 * e + a*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 2] = 8 * e + a*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 3] = 8 * e + a*nns[1] + (b + 1)*ns + c;
					sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
					sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
					sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				}
			}
		}
	}
	

	double *spts= NULL;
	double *sresults=NULL;
	int *seles=NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles=NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{		
		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (double*)malloc(sizeof(double) * 24 * n_bzmesh);
		sresults = (double*)malloc(sizeof(double) * 32 * n_bzmesh);
		seles = (int*)malloc(sizeof(int) * 8 * n_bzmesh);

		displs_spts = (int*)malloc(nProcess * sizeof(int));
		displs_sresults = (int*)malloc(nProcess * sizeof(int));
		displs_seles = (int*)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i<nProcess; i++) {
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 24;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 32;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 8;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 24;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 32;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 8;
		}
	}	

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 8 * 3, MPI_DOUBLE, spts, recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 8 * 4, MPI_DOUBLE, sresults, recvcounts_sresults, displs_sresults ,MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 8, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);
	
	
	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				array<double, 3> pt = { spts[i * 24 + j * 3 + 0], spts[i * 24 + j * 3 + 1], spts[i * 24 + j * 3 + 2] };
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 32 + j * 4 + 0]);
				sresult_all.push_back(sresults[i * 32 + j * 4 + 1]);
				sresult_all.push_back(sresults[i * 32 + j * 4 + 2]);
				sresult_all.push_back(sresults[i * 32 + j * 4 + 3]);
			}
			
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 8> el;
				el[0] = pstart + seles[8 * sum_ele + 0];
				el[1] = pstart + seles[8 * sum_ele + 1];
				el[2] = pstart + seles[8 * sum_ele + 2];
				el[3] = pstart + seles[8 * sum_ele + 3];
				el[4] = pstart + seles[8 * sum_ele + 4];
				el[5] = pstart + seles[8 * sum_ele + 5];
				el[6] = pstart + seles[8 * sum_ele + 6];
				el[7] = pstart + seles[8 * sum_ele + 7];
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 8;
		}
		cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
	}
	
}

void TransportOpt2D::WriteVTK(const vector<array<double, 3>> spt, const vector<double> sdisp, const vector<array<int,8>> sele, int step, string fn)
{
	string fname = fn + "_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout << "POINT_DATA " << sdisp.size() / 4 << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i<sdisp.size() / 4; i++)
		{
			fout << sdisp[i * 4] << " " << sdisp[i * 4 + 1] << " " << sdisp[i * 4 + 2] << "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size() / 4; i++)
		{
			fout << sqrt(sdisp[i * 4] * sdisp[i * 4] + sdisp[i * 4 + 1] * sdisp[i * 4 + 1] + sdisp[i * 4 + 2] * sdisp[i * 4 + 2]) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size() / 4; i++)
		{
			fout << sdisp[4 * i + 3] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void TransportOpt2D::ResultCal_Bezier(double u, double v, const Element2D& bzel, double pt[3], double result[4], double dudx[3], double& detJ)
{
	double dUdx[dim][dim];
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, dim>> dNdx(bzel.IEN.size());
	vector<array<array<double, dim>, dim>> dN2dx2;
	bzel.Para2Phys(u, v, pt);
	BasisFunction(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx, dNdx,dN2dx2, dUdx, detJ);
	result[0] = 0.; result[1] = 0.; result[2] = 0.; result[3] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)	{
		result[0] += Nx[i] * (Vel[dim* bzel.IEN[i] + 0]);
		result[1] += Nx[i] * (Vel[dim* bzel.IEN[i] + 1]);
		result[2] += Nx[i] * (Vel[dim* bzel.IEN[i] + 2]);
		result[3] += Nx[i] * (Pre[bzel.IEN[i]]);
	}
}




