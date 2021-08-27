#include "NS_2Dsteady.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

NS_2Dsteady::NS_2Dsteady()
{
	comm = MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;
}

void NS_2Dsteady::BodyForce(double x, double y, double z, double &Fx, double &Fy, double &Fz)
{
	Fx = 0;
	Fy = 0;
	Fz = 0;
}

void NS_2Dsteady::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open())
	{
		fin_cmat >> neles;
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i < neles; i++)
		{
			if (find(ele_process.begin(), ele_process.end(), i) != ele_process.end())
			{
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;

				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (int j = 0; j < nfunctions; j++)
				{
					for (int k = 0; k < 16; k++)
					{
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else
			{
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
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	add = 0;
	if (fin_bzpt.is_open())
	{
		fin_bzpt >> npts;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++)
		{
			if (find(ele_process.begin(), ele_process.end(), e) != ele_process.end())
			{
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++)
				{
					fin_bzpt >> bzmesh_process[add].pts[i][0] >>
						bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
				}
				add++;
			}
			else
			{
				for (int i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}
		}
		fin_bzpt.close();
		PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

void NS_2Dsteady::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
	case 2:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.2113248654051871;
		Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;
		wght[1] = 1.;
		break;
	}
	case 3:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.1127016653792583;
		Gpt[1] = 0.5;
		Gpt[2] = 0.8872983346207417;
		wght[0] = 0.5555555555555556;
		wght[1] = 0.8888888888888889;
		wght[2] = 0.5555555555555556;
		break;
	}
	case 4:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.06943184420297371;
		Gpt[1] = 0.33000947820757187;
		Gpt[2] = 0.6699905217924281;
		Gpt[3] = 0.9305681557970262;
		wght[0] = 0.3478548451374539;
		wght[1] = 0.6521451548625461;
		wght[2] = 0.6521451548625461;
		wght[3] = 0.3478548451374539;
		break;
	}
	case 5:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.046910077030668;
		Gpt[1] = 0.2307653449471585;
		Gpt[2] = 0.5;
		Gpt[3] = 0.7692346550528415;
		Gpt[4] = 0.953089922969332;
		wght[0] = 0.2369268850561891;
		wght[1] = 0.4786286704993665;
		wght[2] = 0.5688888888888889;
		wght[3] = 0.4786286704993665;
		wght[4] = 0.2369268850561891;
		break;
	}
	default:
	{
		Gpt.resize(2);
		wght.resize(2);
		Gpt[0] = 0.2113248654051871;
		Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;
		wght[1] = 1.;
		break;
	}
	}
}

void NS_2Dsteady::InitializeProblem(const int ndof, const int n_bz, const vector<double> &Vel0, const vector<double> &Pre0, const vector<double> &var)
{

	MPI_Barrier(comm);

	PetscPrintf(PETSC_COMM_WORLD, "Initializing...\n");

	/*Initialize parameters*/
	GaussInfo(4);
	n_bzmesh = n_bz;
	Vel.resize(ndof * dim);
	Pre.resize(ndof);
	Vel = Vel0;
	Pre = Pre0;
	if (var.size() != 0)
	{
		nu = 0.1;
		rou = 0.5;
		alphaM = 0.5 * (3 - rou) / (1 + rou);
		alphaF = 1 / (1 + rou);
		Gama = 0.5 + alphaM - alphaF;
		velocity_max = var[1];
	}
	else
	{
		cerr << "0 variables!\n";
		getchar();
	}
	fx = 0;
	fy = 0;
	fz = 0;

	/*Initialize petsc vector, matrix*/
	PetscInt mat_dim = ndof * dof_all;
	ierr = MatCreate(PETSC_COMM_WORLD, &GK);
	ierr = MatSetSizes(GK, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK, pow(5, dim) * dof_all, NULL, pow(5, dim) * dof_all, NULL);
	MatGetOwnershipRange(GK, &rstart, &rend);
	ierr = MatSetOption(GK, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_solution);
}

void NS_2Dsteady::BasisFunction(double u, double v, int nen, const vector<array<double, 3>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], double &detJ)
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

	double dxdt[dim][dim] = {{0}};
	for (loc = 0; loc < bzpt_num; loc++)
		for (a = 0; a < dim; a++)
			for (b = 0; b < dim; b++)
				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

	double dtdx[dim][dim] = {{0}};
	Matrix2DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i < bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
	}
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.25 * detJ;

	//2nd derivatives
	double dx2dt2[dim][4] = {{0}};
	double dt2dx2[dim][4] = {{0}};
	for (int l = 0; l < dim; l++)
		for (loc = 0; loc < bzpt_num; loc++)
			for (a = 0; a < dim; a++)
				for (b = 0; b < dim; b++)
					dx2dt2[l][dim * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			for (k = 0; k < dim; k++)
				for (a = 0; a < dim; a++)
					for (b = 0; b < dim; b++)
						for (c = 0; c < dim; c++)
							dt2dx2[c][dim * i + j] -= dx2dt2[k][dim * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

	for (loc = 0; loc < bzpt_num; loc++)
		for (i = 0; i < dim; i++)
			for (j = 0; j < dim; j++)
				dN2dx2_bz[loc][i][j] = 0.;

	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (i = 0; i < dim; i++)
		{
			for (j = 0; j < dim; j++)
			{
				for (a = 0; a < dim; a++)
				{
					for (b = 0; b < dim; b++)
					{
						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
					}
					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][dim * i + j];
				}
			}
		}
	}

	for (i = 0; i < nen; i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < dim; m++)
			{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
				for (int n = 0; n < dim; n++)
				{
					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
				}
			}
		}
	}
}

void NS_2Dsteady::PointFormValue(vector<double> &Nx, const vector<array<double, dof_all>> &U, double Value[dof_all])
{
	for (int i = 0; i < dof_all; i++)
		Value[i] = 0;

	for (int i = 0; i < dof_all; i++)
		for (int j = 0; j < Nx.size(); j++)
			Value[i] += U[j][i] * Nx[j];
}

void NS_2Dsteady::PointFormGrad(vector<array<double, dim>> &dNdx, const vector<array<double, dof_all>> &U, double Value[dof_all][dim])
{
	for (int i = 0; i < dof_all; i++)
		for (int j = 0; j < dim; j++)
			Value[i][j] = 0.;

	for (int i = 0; i < dof_all; i++)
		for (int j = 0; j < dim; j++)
			for (int k = 0; k < dNdx.size(); k++)
				Value[i][j] += U[k][i] * dNdx[k][j];
}

void NS_2Dsteady::PointFormHess(vector<array<array<double, dim>, dim>> &d2Ndx2, const vector<array<double, dof_all>> &U, double Value[dof_all][dim][dim])
{
	for (int i = 0; i < dof_all; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				Value[i][j][k] = 0.;
			}
		}
	}
	for (int i = 0; i < dof_all; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				for (int l = 0; l < d2Ndx2.size(); l++)
				{
					Value[i][j][k] += U[l][i] * d2Ndx2[l][j][k];
				}
			}
		}
	}
}

void NS_2Dsteady::Tau(double J[dim][dim], double u[dof_all], double &tauM, double &tauC)
{
	/*calculate stabilization parameter*/
	double C_I = 1.0 / 12.0;

	int i, j, k;

	double G[dim][dim] = {{0}};
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			for (k = 0; k < dim; k++)
				G[i][j] += J[k][i] * J[k][j];

	double g[dim] = {0};
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			g[i] += J[j][i];

	double G_G = 0;
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			G_G += G[i][j] * G[i][j];

	double g_g = 0;
	for (i = 0; i < dim; i++)
		g_g += g[i] * g[i];

	double u_G_u = 0;
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			u_G_u += u[i] * G[i][j] * u[j];

	tauM = u_G_u + C_I * nu * nu * G_G;
	tauM = 1 / sqrt(tauM);

	tauC = (tauM)*g_g;
	tauC = 1 / (tauC);
}

void NS_2Dsteady::FineScale(double tauM, double tauC, double u[dim], double u_x[dim], double u_y[dim], double u_xx[dim], double u_yy[dim], double p, double p_x, double p_y, double u_s[dim], double &p_s)
{
	/*calculate fine scale in VMS*/
	u_s[0] = (u[0] * u_x[0] + u[1] * u_y[0]) + p_x - nu * (u_xx[0] + u_yy[0]) - fx;
	u_s[1] = (u[0] * u_x[1] + u[1] * u_y[1]) + p_y - nu * (u_xx[1] + u_yy[1]) - fy;

	p_s = u_x[0] + u_y[1];

	u_s[0] *= -tauM;
	u_s[1] *= -tauM;
	p_s *= -tauC;
}

void NS_2Dsteady::Residual(vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], const double detJ, const vector<array<double, dof_all>> &U, vector<array<double, dof_all>> Re)
{
	/*calculate residual of the equation*/
	double U_t[dof_all], UU[dof_all];
	double grad_U[dof_all][dim];
	double der2_U[dof_all][dim][dim];

	PointFormValue(Nx, U, UU);
	PointFormGrad(dNdx, U, grad_U);
	PointFormHess(dN2dx2, U, der2_U);

	double u[dim] = {UU[0], UU[1]};
	double u_x[dim] = {grad_U[0][0], grad_U[1][0]};
	double u_y[dim] = {grad_U[0][1], grad_U[1][1]};

	double u_xx[dim] = {der2_U[0][0][0], der2_U[1][0][0]};
	double u_yy[dim] = {der2_U[0][1][1], der2_U[1][1][1]};

	double p = UU[dim];
	double p_x = grad_U[dim][0], p_y = grad_U[dim][1];

	double InvGradMap[dim][dim];
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			InvGradMap[i][j] = dudx[i][j];
		}
	}

	double tauM, tauC;
	Tau(InvGradMap, UU, tauM, tauC);
	double u_s[dim], p_s;
	FineScale(tauM, tauC, u, u_x, u_y, u_xx, u_yy, p, p_x, p_y, u_s, p_s);

	int a, nen = Nx.size();
	for (a = 0; a < nen; a++)
	{
		double Na = Nx[a];
		double Na_x = dNdx[a][0];
		double Na_y = dNdx[a][1];

		double Rux, Ruy, Rp;

		Rux = -Na * fx;
		Ruy = -Na * fy;
		Rp = 0.0;

		Rux += -Na_x * p + nu * (Na_x * (u_x[0] + u_x[0]) + Na_y * (u_y[0] + u_x[1]));
		Ruy += -Na_y * p + nu * (Na_x * (u_x[1] + u_y[0]) + Na_y * (u_y[1] + u_y[1]));
		Rp += Na * (u_x[0] + u_y[1]);

		Rux += -(Na_x * p_s);
		Ruy += -(Na_y * p_s);
		Rp += -(Na_x * u_s[0] + Na_y * u_s[1]);

		Rux += +Na * ((u[0] + u_s[0]) * u_x[0] + (u[1] + u_s[1]) * u_y[0]);
		Ruy += +Na * ((u[0] + u_s[0]) * u_x[1] + (u[1] + u_s[1]) * u_y[1]);

		Rux += -(Na_x * u_s[0] * (u[0] + u_s[0]) + Na_y * u_s[0] * (u[1] + u_s[1]));
		Ruy += -(Na_x * u_s[1] * (u[0] + u_s[0]) + Na_y * u_s[1] * (u[1] + u_s[1]));

		Re[a][0] += Rux * detJ;
		Re[a][1] += Ruy * detJ;
		Re[a][2] += Rp * detJ;
	}
}

void NS_2Dsteady::Tangent(vector<double> &Nx, vector<array<double, dim>> &dNdx, double dudx[dim][dim], const double detJ, const vector<array<double, dof_all>> &U, vector<array<vector<array<double, dof_all>>, dof_all>> &Ke)
{
	/*calculate tangent matrix*/
	double u[dof_all];
	PointFormValue(Nx, U, u);
	double ux = u[0];
	double uy = u[1];

	double InvGradMap[dim][dim];
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			InvGradMap[i][j] = dudx[i][j];
		}
	}

	double tauM, tauC;
	Tau(InvGradMap, u, tauM, tauC);
	int a, b, nen = Nx.size();
	for (a = 0; a < nen; a++)
	{
		double Na = Nx[a];
		double Na_x = dNdx[a][0];
		double Na_y = dNdx[a][1];
		for (b = 0; b < nen; b++)
		{
			double Nb = Nx[b];
			double Nb_x = dNdx[b][0];
			double Nb_y = dNdx[b][1];
			/* ----- */
			int i, j;
			double T[dof_all][dof_all];
			double Tii =
				(+Na * (ux * Nb_x + uy * Nb_y) + nu * (Na_x * Nb_x + Na_y * Nb_y) + tauM * (ux * Na_x + uy * Na_y) *
																						/**/ ((ux * Nb_x + uy * Nb_y)));
			T[0][0] = (+nu * Na_x * Nb_x + tauC * Na_x * Nb_x);
			T[0][1] = (+nu * Na_y * Nb_x + tauC * Na_x * Nb_y);
			//
			T[1][0] = (+nu * Na_x * Nb_y + tauC * Na_y * Nb_x);
			T[1][1] = (+nu * Na_y * Nb_y + tauC * Na_y * Nb_y);
			//
			T[0][0] += Tii;
			T[1][1] += Tii;

			T[0][2] = -Na_x * Nb + tauM * (ux * Na_x + uy * Na_y) * Nb_x;
			T[1][2] = -Na_y * Nb + tauM * (ux * Na_x + uy * Na_y) * Nb_y;

			T[2][0] = +Na * Nb_x + tauM * Na_x * ((ux * Nb_x + uy * Nb_y));
			T[2][1] = +Na * Nb_y + tauM * Na_y * ((ux * Nb_x + uy * Nb_y));

			T[2][2] = +tauM * (Na_x * Nb_x + Na_y * Nb_y);

			for (i = 0; i < dof_all; i++)
				for (j = 0; j < dof_all; j++)
					Ke[a][i][b][j] += T[i][j] * detJ;
		}
	}
}

void NS_2Dsteady::MatrixAssembly(vector<array<vector<array<double, dof_all>>, dof_all>> &Ke, const vector<int> &IEN, Mat &GK)
{
	int i, j, A, B, m, n;
	int row_start, row_end, row_now;
	int add = 0;

	PetscInt *nodeList = new PetscInt[IEN.size() * dof_all];
	PetscReal *tmpGK = new PetscReal[IEN.size() * dof_all * IEN.size() * dof_all];

	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		for (i = 0; i < dof_all; i++)
		{
			nodeList[dof_all * m + i] = dof_all * A + i;
			for (n = 0; n < IEN.size(); n++)
			{
				B = IEN[n];
				for (j = 0; j < dof_all; j++)
				{
					tmpGK[add] = Ke[m][i][n][j];
					add++;
				}
			}
		}
	}
	MatSetValues(GK, IEN.size() * dof_all, nodeList, IEN.size() * dof_all, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;
}

void NS_2Dsteady::ResidualAssembly(vector<array<double, dof_all>> &Re, const vector<int> &IEN, Vec &GR)
{
	int i, j, A, B, m, n;
	int add = 0;

	PetscInt *nodeList = new PetscInt[IEN.size() * dof_all];
	PetscReal *tmpGR = new PetscReal[IEN.size() * dof_all];
	for (m = 0; m < IEN.size(); m++)
	{
		A = IEN[m];
		for (i = 0; i < dof_all; i++)
		{
			nodeList[dof_all * m + i] = dof_all * A + i;
			tmpGR[add] = -Re[m][i];
			add++;
		}
	}

	VecSetValues(GR, IEN.size() * dof_all, nodeList, tmpGR, ADD_VALUES);
	delete nodeList;
	delete tmpGR;
}

void NS_2Dsteady::BuildLinearSystemProcess(const vector<Vertex2D> &cpts, const vector<array<double, dim>> &velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node)
{
	/*Build linear system in each process*/
	int e;
	cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for " << bzmesh_process.size() << " elements.\n";
	for (e = 0; e < bzmesh_process.size(); e++)
	{
		int nen, A;
		double ux_bc, uy_bc, p_bc;

		double dudx[dim][dim];
		double detJ;
		vector<double> Nx;
		vector<array<double, dim>> dNdx;
		vector<array<array<double, dim>, dim>> dN2dx2;
		vector<array<double, dof_all>> Re;
		vector<array<vector<array<double, dof_all>>, dof_all>> Ke;

		vector<array<double, dof_all>> v_node;

		nen = bzmesh_process[e].IEN.size();

		Nx.clear();
		Nx.resize(nen, 0);
		dNdx.clear();
		dNdx.resize(nen, {0});
		dN2dx2.clear();
		dN2dx2.resize(nen, {{0}});
		Re.clear();
		Re.resize(nen);
		Ke.clear();
		Ke.resize(nen);
		v_node.clear();
		v_node.resize(nen);

		for (int i = 0; i < nen; i++)
			for (int j = 0; j < dof_all; j++)
				Ke[i][j].resize(nen);

		for (int i = 0; i < nen; i++)
		{
			for (int m = 0; m < dof_all; m++)
			{
				for (int j = 0; j < nen; j++)
				{
					for (int n = 0; n < dof_all; n++)
					{
						Ke[i][m][j][n] = 0.;
					}
				}
				Re[i][m] = 0.;
			}
		}

		for (int i = 0; i < nen; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				v_node[i][j] = velocity_max * velocity_node[bzmesh_process[e].IEN[i] * dim + j];
			}
			v_node[i][dim] = pressure_node[bzmesh_process[e].IEN[i]];
		}
		// cout << "Integrating...\n";
		for (int i = 0; i < Gpt.size(); i++)
		{
			for (int j = 0; j < Gpt.size(); j++)
			{
				BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				Tangent(Nx, dNdx, dudx, detJ, v_node, Ke);
				Residual(Nx, dNdx, dN2dx2, dudx, detJ, v_node, Re);
			}
		}

		/*Apply Boundary Condition*/
		for (int i = 0; i < nen; i++)
		{
			A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1)
			{
				//inlet
				ux_bc = velocity_max * velocity_bc[A][0] - velocity_node[A * dim];
				ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
				uy_bc = velocity_max * velocity_bc[A][1] - velocity_node[A * dim + 1];
				ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
				// cout << "ux_bc = " << ux_bc << 	"uy_bc = " << uy_bc	<< endl;
			}
			else if (cpts[A].label == 0)
			{
				//wall
				ux_bc = 0.0 - velocity_node[A * dim];
				ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
				uy_bc = 0.0 - velocity_node[A * dim + 1];
				ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
			}
			// else if (cpts[A].label == 4)
			// {
			// 	wall
			// 	ux_bc = -5 * velocity_max * velocity_bc[A][0] - velocity_node[A * dim];
			// 	ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
			// 	uy_bc = -5 * velocity_max * velocity_bc[A][1] - velocity_node[A * dim + 1];
			// 	ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
			// }
			// else if (cpts[A].label == 2 || cpts[A].label == 3)
			// {
			// 	ux_bc = velocity_max * 0.12 * velocity_bc[A][0] - velocity_node[A * dim];
			// 	ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
			// 	uy_bc = velocity_max * 0.12 * velocity_bc[A][1] - velocity_node[A * dim + 1];

			// 	ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
			// 	// p_bc = 0.0 - pressure_node[A];
			// 	// ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label >= 4 && cpts[A].label <= 9)
			// {
			// 	ux_bc = velocity_max * 0.08 * velocity_bc[A][0] - velocity_node[A * dim];
			// 	ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
			// 	uy_bc = velocity_max * 0.08 * velocity_bc[A][1] - velocity_node[A * dim + 1];
			// 	ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
			// }
			// else if (cpts[A].label == 11)
			// {
			// 	ux_bc = velocity_max * 0.05 * velocity_bc[A][0] - velocity_node[A * dim];
			// 	ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
			// 	uy_bc = velocity_max * 0.05 * velocity_bc[A][1] - velocity_node[A * dim + 1];
			// 	ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
			// }
			// if (cpts[A].label >= 2)
			// {
			// 	p_bc = 0.0 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }

			// if (cpts[A].label == 2)
			// {
			// 	p_bc = -0.65 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 3)
			// {
			// 	p_bc = -.65 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 4)
			// {
			// 	p_bc = -1.65 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 5)
			// {
			// 	p_bc = -1.65 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 6)
			// {
			// 	p_bc = -1.65 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 7)
			// {
			// 	p_bc = -1.25 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 8)
			// {
			// 	p_bc = -1.2 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 9)
			// {
			// 	p_bc = -1.2 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 10)
			// {
			// 	p_bc = -0.25 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
			// else if (cpts[A].label == 11)
			// {
			// 	p_bc = -0.15 - pressure_node[A];
			// 	ApplyBoundaryCondition(p_bc, i, 2, Ke, Re);
			// }
		}

		/*Start element vector assembly*/
		ResidualAssembly(Re, bzmesh_process[e].IEN, GR);

		/*Start element matrix assembly*/
		MatrixAssembly(Ke, bzmesh_process[e].IEN, GK);
	}
	cout << "Process " << comRank << " :complete build matrix and vector!\n";
	VecAssemblyBegin(GR);
	MatAssemblyBegin(GK, MAT_FINAL_ASSEMBLY);
}

void NS_2Dsteady::ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, dof_all>>, dof_all>> &Ke, vector<array<double, dof_all>> &Re)
{
	int j, k;
	for (j = 0; j < Re.size(); j++)
	{
		for (k = 0; k < dof_all; k++)
		{
			Re[j][k] += bc_value * Ke[j][k][pt_num][variable_num];
		}
	}
	for (j = 0; j < Re.size(); j++)
	{
		for (k = 0; k < dof_all; k++)
		{
			Ke[pt_num][variable_num][j][k] = 0.0;
			Ke[j][k][pt_num][variable_num] = 0.0;
		}
	}
	Re[pt_num][variable_num] = -bc_value;
	Ke[pt_num][variable_num][pt_num][variable_num] = 1.0;
}

void NS_2Dsteady::AssignProcessor(vector<vector<int>> &ele_proc)
{
	/*Assign the partitioned bezier elements to this processor */
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void NS_2Dsteady::Run(const vector<Vertex2D> &cpts, const vector<Element2D> &tmesh, const vector<array<double, dim>> &velocity_bc, string fn)
{

	int n_iterator(3);
	int l(0);
	time_t t0, t1;
	vector<double> V_delta(dim * cpts.size()), P_delta(cpts.size());

	ReadBezierElementProcess(fn);
	for (l = 0; l < n_iterator; l++)
	{

		PetscPrintf(PETSC_COMM_WORLD, "Iteration Step: %d\n", l);
		/*Build Linear System*/
		PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");
		MatZeroEntries(GK);
		VecSet(GR, 0);
		t0 = time(NULL);
		BuildLinearSystemProcess(cpts, velocity_bc, Vel, Pre);
		VecAssemblyEnd(GR);
		PetscPrintf(PETSC_COMM_WORLD, "Done Vector Assembly...\n");
		MatAssemblyEnd(GK, MAT_FINAL_ASSEMBLY);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 - t0);
		MPI_Barrier(comm);
		PetscPrintf(PETSC_COMM_WORLD, "Solving...\n");

		DebugVisualizeVec(GR, fn + "debug/GR" + to_string(l) + ".m");
		DebugVisualizeMat(GK, fn + "debug/GK" + to_string(l) + ".m");

		t0 = time(NULL);
		/*Petsc KSP solver setting*/
		KSPCreate(PETSC_COMM_WORLD, &ksp);
		KSPSetOperators(ksp, GK, GK);
		KSPSetFromOptions(ksp);

		// KSPGetPC(ksp, &pc);
		// PCSetType(pc, PCBJACOBI);
		// KSPSetType(ksp, KSPGMRES);
		// KSPGMRESSetRestart(ksp, 500);
		// KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
		// KSP *subksp;
		// PC subpc;
		// PetscInt first,nlocal;
		// KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, 100000);
		// KSPSetPCSide(ksp, PC_RIGHT);

		// KSPSetUp(ksp);
		// PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
		// for (int i = 0; i<nlocal; i++) {
		// 	KSPGetPC(subksp[i], &subpc);
		// 	PCSetType(subpc, PCILU);
		// 	KSPSetType(subksp[i], KSPGMRES);
		// 	KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE);
		// 	KSPSetPCSide(subksp[i], PC_RIGHT);
		// }

		// string work_dir = "../io/single_pipe4/";
		// DebugVisualizeMat(GK, work_dir + "debug/NS_GK_afterBC.m");
		// DebugVisualizeVec(GR, work_dir + "debug/NS_GR_afterBC.m");

		// for (uint i = 0; i < cpts.size(); i++){
		// 	cout << Vel[dim * i] << " " << Vel[dim * i + 1 ] << " " << 0 << "\n";
		// }

		/*Solving the equation*/
		KSPSolve(ksp, GR, temp_solution);
		// KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "iterations %d\n", its);
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
		PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Solving with time: %d \n", t1 - t0);

		PetscScalar GRVal, sol_Val;
		VecNorm(GR, NORM_2, &GRVal);
		VecNorm(temp_solution, NORM_2, &sol_Val);
		PetscPrintf(PETSC_COMM_WORLD, "Residual: %lf \n", GRVal);
		PetscPrintf(PETSC_COMM_WORLD, "|dx| = %f \n", sol_Val);

		if (sol_Val < 1e-8 || GRVal < 1e-8)
		{
			break;
		}

		/*Collect the solution from all processors*/
		Vec temp_solution_seq;
		VecScatter ctx;
		PetscScalar *_a;
		VecScatterCreateToAll(temp_solution, &ctx, &temp_solution_seq);
		VecScatterBegin(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecGetArray(temp_solution_seq, &_a);
		MPI_Barrier(comm);

		for (uint i = 0; i < cpts.size(); i++)
		{
			V_delta[dim * i] = PetscRealPart(_a[dof_all * i]);
			V_delta[dim * i + 1] = PetscRealPart(_a[dof_all * i + 1]);
			P_delta[i] = PetscRealPart(_a[dof_all * i + dim]);
		}
		VecRestoreArray(temp_solution_seq, &_a);
		VecScatterDestroy(&ctx);
		VecDestroy(&temp_solution_seq);

		for (uint i = 0; i < cpts.size(); i++)
		{
			Vel[dim * i] += V_delta[dim * i];
			Vel[dim * i + 1] += V_delta[dim * i + 1];
			Pre[i] += P_delta[i];
		}
		MPI_Barrier(comm);
		/*Visualize the result*/
		if (comRank == 0)
		{
			cout << "Visualizing...\n";
			VisualizeVTK_ControlMesh(cpts, tmesh, l, fn);
		}
		MPI_Barrier(comm);
		VisualizeVTK_PhysicalDomain(l, fn + "final_physics");
	}
	MatDestroy(&GK);
	VecDestroy(&GR);
	VecDestroy(&temp_solution);
	KSPDestroy(&ksp);

	MPI_Barrier(comm);
}

void NS_2Dsteady::VisualizeVTK_ControlMesh(const vector<Vertex2D> &spt, const vector<Element2D> &mesh, int step, string fn)
{
	stringstream ss;
	ss << step;

	string fname = fn + "controlmesh_VelocityPressure_" + ss.str() + ".vtk";
	// string fname = fn + "controlmesh_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << Vel.size() / dim << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i < Vel.size() / dim; i++)
		{
			fout << Vel[i * dim] << " " << Vel[i * dim + 1] << " " << 0.0 << "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < Pre.size(); i++)
		{
			fout << sqrt(Vel[i * dim] * Vel[i * dim] + Vel[i * dim + 1] * Vel[i * dim + 1]) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < Pre.size(); i++)
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
		for (uint i = 0; i < Vel.size() / dim; i++)
		{
			fout1 << Vel[i * dim] << " " << Vel[i * dim + 1] << " " << 0.0 << "\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void NS_2Dsteady::VisualizeVTK_PhysicalDomain(int step, string fn)
{
	vector<array<double, 3>> spt_all; // sample points
	vector<double> sresult_all;
	vector<array<int, 4>> sele_all;

	const int ns_ele = 4;
	const int ns_pt = 9;
	double detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	double spt_proc[num_bzmesh_ele * ns_pt * 3];
	double sresult_proc[num_bzmesh_ele * ns_pt * dof_all];
	int sele_proc[num_bzmesh_ele * ns_ele * 4];
	for (unsigned int e = 0; e < num_bzmesh_ele; e++)
	{
		int ns(3);
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a < ns; a++)
		{
			for (int b = 0; b < ns; b++)
			{
				double pt1[phy_dim], dudx[dim];
				double result[dof_all];
				ResultCal_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
				spt_proc[e * ns_pt * 3 + loc * 3 + 0] = pt1[0];
				spt_proc[e * ns_pt * 3 + loc * 3 + 1] = pt1[1];
				spt_proc[e * ns_pt * 3 + loc * 3 + 2] = pt1[2];
				for (int c = 0; c < dof_all; c++)
				{
					sresult_proc[dof_all * ns_pt * e + loc * dof_all + c] =
						result[c];
				}
				loc++;
			}
		}
		int nns[2] = {ns * ns, ns};
		loc = 0;
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				sele_proc[e * ns_ele * 4 + loc * 4 + 0] = e * ns_pt + a * ns + b;
				sele_proc[e * ns_ele * 4 + loc * 4 + 1] = e * ns_pt + a * ns + b + 1;
				sele_proc[e * ns_ele * 4 + loc * 4 + 2] =
					e * ns_pt + (a + 1) * ns + (b + 1);
				sele_proc[e * ns_ele * 4 + loc * 4 + 3] = e * ns_pt + (a + 1) * ns + b;
				// cout << "e: " << e << " loc: " << loc << " " << sele_proc[e * ns_ele
				// * 4 + loc * 4 + 0] << " " << sele_proc[e * ns_ele * 4 + loc * 4 + 1]
				// << " " << sele_proc[e * ns_ele * 4 + loc * 4 + 2] << " " <<
				// sele_proc[e * ns_ele * 4 + loc * 4 + 3] << endl;
				loc++;
			}
		}
	}

	double *spts = NULL;
	double *sresults = NULL;
	int *seles = NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles = NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{
		num_bzmesh_eles = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_spts = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_sresults = (int *)malloc(sizeof(int) * nProcess);
		recvcounts_seles = (int *)malloc(sizeof(int) * nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (double *)malloc(sizeof(double) * 3 * n_bzmesh * ns_pt);
		sresults = (double *)malloc(sizeof(double) * dof_all * n_bzmesh * ns_pt);
		seles = (int *)malloc(sizeof(int) * 4 * n_bzmesh * ns_ele);

		displs_spts = (int *)malloc(nProcess * sizeof(int));
		displs_sresults = (int *)malloc(nProcess * sizeof(int));
		displs_seles = (int *)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i < nProcess; i++)
		{
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 3 * ns_pt;
			displs_sresults[i] =
				displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * dof_all * ns_pt;
			displs_seles[i] =
				displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4 * ns_ele;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * ns_pt * 3;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * ns_pt * dof_all;
			recvcounts_seles[i] = num_bzmesh_eles[i] * ns_ele * 4;
		}
	}

	MPI_Gatherv(spt_proc, num_bzmesh_ele * ns_pt * 3, MPI_DOUBLE, spts, recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * ns_pt * dof_all, MPI_DOUBLE, sresults, recvcounts_sresults, displs_sresults, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * ns_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < ns_pt; j++)
			{
				array<double, 3> pt = {spts[i * ns_pt * 3 + j * 3 + 0],
									   spts[i * ns_pt * 3 + j * 3 + 1],
									   spts[i * ns_pt * 3 + j * 3 + 2]};
				spt_all.push_back(pt);
				for (int c = 0; c < dof_all; c++)
				{
					sresult_all.push_back(
						sresults[i * dof_all * ns_pt + j * dof_all + c]);
				}
			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				for (int nse = 0; nse < ns_ele; nse++)
				{
					array<int, 4> el;
					el[0] = pstart + seles[ns_ele * sum_ele * 4 + nse * 4 + 0];
					el[1] = pstart + seles[ns_ele * sum_ele * 4 + nse * 4 + 1];
					el[2] = pstart + seles[ns_ele * sum_ele * 4 + nse * 4 + 2];
					el[3] = pstart + seles[ns_ele * sum_ele * 4 + nse * 4 + 3];
					sele_all.push_back(el);
				}
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * ns_pt;
		}
		// cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
	}

	if (spts)
		free(spts);
	if (sresults)
		free(sresults);
	if (seles)
		free(seles);
	if (displs_spts)
		free(displs_spts);
	if (displs_sresults)
		free(displs_sresults);
	if (displs_seles)
		free(displs_seles);
	if (num_bzmesh_eles)
		free(num_bzmesh_eles);
	if (recvcounts_spts)
		free(recvcounts_spts);
	if (recvcounts_sresults)
		free(recvcounts_sresults);
	if (recvcounts_seles)
		free(recvcounts_seles);
}

void NS_2Dsteady::WriteVTK(const vector<array<double, 3>> spt, const vector<double> sdisp, const vector<array<int, 4>> sele, int step, string fn)
{
	string fname = fn + "_VelocityPressure" + "_" + to_string(step) + ".vtk";
	// string fname = fn + "_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i < spt.size(); i++)
		{
			fout << std::setprecision(9) << spt[i][0] << std::fixed << " "
				 << std::setprecision(9) << spt[i][1] << std::fixed << " "
				 << std::setprecision(9) << spt[i][2] << std::fixed << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i < sele.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << sdisp.size() / dof_all << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i < sdisp.size() / dof_all; i++)
		{
			fout << sdisp[i * dof_all] << " " << sdisp[i * dof_all + 1] << " " << 0.0 << "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sdisp.size() / dof_all; i++)
		{
			fout << sqrt(sdisp[i * dof_all] * sdisp[i * dof_all] + sdisp[i * dof_all + 1] * sdisp[i * dof_all + 1]) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i < sdisp.size() / dof_all; i++)
		{
			fout << sdisp[dof_all * i + 2] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NS_2Dsteady::ResultCal_Bezier(double u, double v, const Element2D &bzel, double pt[phy_dim], double result[dof_all], double dudx[dim], double &detJ)
{
	double dUdx[dim][dim];
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, dim>> dNdx(bzel.IEN.size());
	vector<array<array<double, dim>, dim>> dN2dx2;
	bzel.Para2Phys(u, v, pt);
	BasisFunction(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx, dNdx, dN2dx2, dUdx, detJ);
	result[0] = 0.;
	result[1] = 0.;
	result[2] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		result[0] += Nx[i] * (Vel[dim * bzel.IEN[i] + 0]);
		result[1] += Nx[i] * (Vel[dim * bzel.IEN[i] + 1]);
		result[2] += Nx[i] * (Pre[bzel.IEN[i]]);
	}
}
