#include "LBO.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

LBO::LBO()
{
}

void LBO::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch(ng)
	{
	case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;							wght[1]=1.;
			break;
		}
	case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.1127016653792583;			Gpt[1]=0.5;							Gpt[2]=0.8872983346207417;
			wght[0]=0.5555555555555556;			wght[1]=0.8888888888888889;			wght[2]=0.5555555555555556;
			break;
		}
	case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.06943184420297371;			Gpt[1]=0.33000947820757187;			Gpt[2]=0.6699905217924281;			Gpt[3]=0.9305681557970262;
			wght[0]=0.3478548451374539;			wght[1]=0.6521451548625461;			wght[2]=0.6521451548625461;			wght[3]=0.3478548451374539;
			break;
		}
	default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;							wght[1]=1.;
			break;
		}
	}
}

void LBO::InitializeProblem(const vector<double>& N_plus0, const vector<double>& N_minus0, const vector<double>& var, double tstep)
{
	//UpdateCA(CA0);
	N_plus.resize(N_plus0.size());
	N_minus.resize(N_minus0.size());
	par.resize(var.size());
	N_plus = N_plus0;
	N_minus = N_minus0;
	par = var;
	dt = tstep;
}

void LBO::UpdateCA(vector<double>& CA0)
{
		CA_B.clear();
		CA_B.resize(CA0.size());
		CA_B = CA0;
}

void LBO::BasisFunction(double u, double v, const vector<array<double, 3>>& pt, double Nx[4], double dNdx[4][3], double& detJ)
{
	double Nu[2]={(1.-u), u};
	double Nv[2]={(1.-v), v};
	double dNdu[2]={-1., 1.};
	double dNdv[2] = { -1., 1. };
	double dNdt[4][2];
	int i,j,a,b,loc(0);
	for (j = 0; j < 2; j++)
	{
		Nx[loc] = Nu[j] * Nv[0];
		dNdt[loc][0] = dNdu[j] * Nv[0];
		dNdt[loc][1] = Nu[j] * dNdv[0];
		loc++;
	}
	for (j= 1; j >= 0; j--)
	{
		Nx[loc] = Nu[j] * Nv[1];
		dNdt[loc][0] = dNdu[j] * Nv[1];
		dNdt[loc][1] = Nu[j] * dNdv[1];
		loc++;
	}

	double dxdt[3][2]={{0.,0.},{0.,0.},{0.,0.}};
	loc=0;
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			for (a=0;a<3;a++)
			{
				for(b=0;b<2;b++)
				{
					dxdt[a][b]+=pt[loc][a]*dNdt[loc][b];
				}
			}
			loc++;
		}
	}
	double g[2][2]={{0.,0.},{0.,0.}};
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			g[i][j]=dxdt[0][i]*dxdt[0][j]+dxdt[1][i]*dxdt[1][j]+dxdt[2][i]*dxdt[2][j];
		}
	}
	double gdet=g[0][0]*g[1][1]-g[0][1]*g[1][0];
	double gi[2][2]={{g[1][1]/gdet,-g[0][1]/gdet},{-g[1][0]/gdet,g[0][0]/gdet}};
	for(i=0;i<4;i++)
	{
		for(a=0;a<3;a++)
		{
			dNdx[i][a]=(gi[0][0]*dNdt[i][0]+gi[0][1]*dNdt[i][1])*dxdt[a][0]+(gi[1][0]*dNdt[i][0]+gi[1][1]*dNdt[i][1])*dxdt[a][1];
		}
	}
	detJ=0.25*sqrt(g[0][0]*g[1][1]-g[0][1]*g[1][0]); 
}

void LBO::ElementMassMatrix(double Nx[4], double detJ, double EM[4][4])
{
	int i, j;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			EM[i][j] += Nx[i] * Nx[j] * detJ;
		}
	}
}

void LBO::ElementStiffMatrix(double dNdx[4][3], double detJ, double EK[4][4])
{
	int i, j;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void LBO::ElementConvectionMatrix(double Nx[4], double dNdx[4][3], double v[3], double detJ, double EC[4][4])
{
	int i, j;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			EC[i][j] += Nx[i] * (v[0] * dNdx[j][0] + v[1] * dNdx[j][1] + v[2] * dNdx[j][2])*detJ;
		}
	}
}

void LBO::Assembly(double EM[4][4], double EC_plus[4][4], double EC_minus[4][4], const vector<int>& IEN, SparseMatrix<double>& GM, SparseMatrix<double>& GC_plus, SparseMatrix<double>& GC_minus)
{
	unsigned int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			GM.coeffRef(A, B) += EM[i][j];
			GC_plus.coeffRef(A, B) += EC_plus[i][j];
			GC_minus.coeffRef(A, B) += EC_minus[i][j];
		}
	}
}

void LBO::BuildLinearSystem(const vector<Element2D>& mesh, const double Vplus, const double Vminus, SparseMatrix<double>& GM, SparseMatrix<double>& GC_plus, SparseMatrix<double>& GC_minus)
{
	unsigned int e, i, j, k, a, b;
	double EM[4][4], EC_plus[4][4], EC_minus[4][4];
	double Nx[4];
	double dNdx[4][3];
	double detJ, ca;
	
	for (e = 0; e<mesh.size(); e++)
	{
		double v_plus[3] = { 0.,0.,Vplus };
		double v_minus[3] = { 0.,0.,Vminus };
		for (i = 0; i<4; i++)
		{
			for (j = 0; j<4; j++)
			{
				EM[i][j] = 0.;
				EC_plus[i][j] = 0.;
				EC_minus[i][j] = 0.;
			}
		}
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				BasisFunction(Gpt[i], Gpt[j], mesh[e].pts, Nx, dNdx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				ElementMassMatrix(Nx, detJ, EM);
				ElementConvectionMatrix(Nx,dNdx,v_plus, detJ, EC_plus);
				ElementConvectionMatrix(Nx, dNdx, v_minus, detJ, EC_minus);
				
			}
		}
		Assembly(EM, EC_plus, EC_minus, mesh[e].cnct, GM, GC_plus, GC_minus);

		//string fn44, fn45, fn46;
		//fn44 = "../io/matrix_examine/LBO_Element_MassMatrix.txt";
		//fn45 = "../io/matrix_examine/LBO_Element_StiffMatrix.txt";
		//ofstream fout44, fout45, fout46;
		//fout44.open(fn44, ios::app);
		//fout45.open(fn45, ios::app);
		//fout44 << e << "\n";
		//fout45 << e << "\n";
		//for (i = 0; i<4; i++)
		//{
		//	for (j = 0; j<4; j++)
		//	{
		//		fout44 << EM[i][j] << " ";
		//		fout45 << EK[i][j] << " ";
		//
		//	}
		//	fout44 << "\n";
		//	fout45 << "\n";
		//}
		//fout44 << "\n";
		//fout45 << "\n";
		//fout44.close();
		//fout45.close();

	}
}

void LBO::Solver(SparseMatrix<double>& GM, SparseMatrix<double>& GC_plus, SparseMatrix<double>& GC_minus, vector<double>& Bv, vector<double>& CA)
{
	Bv.clear();
	Bv.resize(N_plus.size());
	CA.clear();
	CA.resize(CA_B.size());
	//SparseMatrix<double> GK_CA = GM - dt*par[0] * GK;
	//SparseMatrix<double> GK_X = GM - dt*par[1] * GK;
	//SparseMatrix<double> GK_B = GM - dt*par[2] * GK;
	SparseMatrix<double> GM_plus = GM +dt*GC_plus;
	SparseMatrix<double> GM_minus = GM+dt*GC_minus;
	
	VectorXd tmp_Nplus(N_plus.size()), tmp_Nminus(N_minus.size()),tmp_CA(CA_B.size());
	for (uint i = 0; i < N_plus.size(); i++)
	{
		tmp_CA[i] = CA_B[i];
	}
	for (uint i = 0; i < N_plus.size(); i++)
	{
		tmp_Nplus[i] = N_plus[i];
	}
	for (uint i = 0; i < N_minus.size(); i++)
	{
		tmp_Nminus[i] = N_minus[i];
	}
	VectorXd GF_plus = GM*tmp_Nplus;
	VectorXd GF_minus = GM*tmp_Nminus;
	VectorXd tmp_plus = par[3] * (GM*tmp_CA) - par[5] * (GM*tmp_Nplus);
	VectorXd tmp_minus = par[4] * (GM*tmp_CA) - par[6] * (GM*tmp_Nminus);
	for (uint i = 0; i < Bv.size(); i++)
	{
		Bv[i] = tmp_plus[i]+tmp_minus[i];
		GF_plus[i] += dt*tmp_plus[i];
		GF_minus[i] += dt*tmp_minus[i];
	}

	for (uint i = 0; i < 121; i++)
	{
		GM_plus.coeffRef(i, i) *= 1.e11;
		GM_minus.coeffRef(i, i) *= 1.e11;
		GF_plus[i] = 0.6 * GM_plus.coeffRef(i, i);
		GF_minus[i] = 0 * GM_minus.coeffRef(i, i);
	}
	for (uint i = 481; i < 602; i++)
	{
		GM_plus.coeffRef(i, i) *= 1.e11;
		GM_minus.coeffRef(i, i) *= 1.e11;
		GF_plus[i] = 0.18 * GM_plus.coeffRef(i, i);
		GF_minus[i] = 0 * GM_minus.coeffRef(i, i);
	}
	//string fn4;
	//fn4 = "../io/matrix_examine/LBO_GM.txt";
	//ofstream fout4;
	//fout4.open(fn4);
	//fout4 << GM;
	//fout4.close();
	//
	//string fn5;
	//fn5 = "../io/matrix_examine/LBO_GK.txt";
	//ofstream fout5;
	//fout5.open(fn5);
	//fout5 << GK;
	//fout5.close();
	//
	//string fn6;
	//fn6 = "../io/matrix_examine/LBO_Bv.txt";
	//ofstream fout6;
	//fout6.open(fn6);
	//for (int i = 0; i<N_plus.size(); i++)
	//	fout6 << i <<" "<<Bv[i] << "\n";
	//fout6.close();
	
	//string fn7;
	//fn7 = "../io/matrix_examine/LBO_GFB.txt";
	//ofstream fout7;
	//fout7.open(fn7);
	//fout7 << GF_B;
	//fout7.close();

	//SimplicialLDLT<SparseMatrix<double>> solver1;
	SimplicialLDLT<SparseMatrix<double>> solver2;
	SimplicialLDLT<SparseMatrix<double>> solver3;
	//solver1.compute(GM + dt*par[0] * GK);
	solver2.compute(GM_plus);
	solver3.compute(GM_minus);
	//VectorXd sol_CA = solver1.solve(GF_CA);
	VectorXd sol_plus = solver2.solve(GF_plus);
	VectorXd sol_minus = solver3.solve(GF_minus);
	//for (uint i = 0; i<CA_B.size(); i++)
	//{
	//	CA_B[i] = sol_CA[i];
	//	CA[i] = sol_CA[i];
	//}
	for (uint i = 0; i<N_plus.size(); i++)
	{
		N_plus[i] = sol_plus[i];
	}
	for (uint i = 0; i<N_minus.size(); i++)
	{
		N_minus[i] = sol_minus[i];
	}
}

void LBO::VisualizeVTK(const vector<array<double, 3>>& spt, const vector<Element2D>& mesh, string fn)
{
	string fname = fn + "_N_plus.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].cnct[0] << " " << mesh[i].cnct[1] << " " << mesh[i].cnct[2] << " " << mesh[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << N_plus.size() << "\nSCALARS N_plus float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<N_plus.size(); i++)
		{
			fout << N_plus[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	fname = fn + "_N_minus.vtk";
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].cnct[0] << " " << mesh[i].cnct[1] << " " << mesh[i].cnct[2] << " " << mesh[i].cnct[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << N_minus.size() << "\nSCALARS N_minus float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<N_minus.size(); i++)
		{
			fout << N_minus[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void LBO::Run(const vector<array<double, 3>>& pts, const vector<Element2D>& mesh, vector<double>& CA, vector<double>& Bv, string fn)
{
	cout << "Reaction...\n";
	GaussInfo(3);
	UpdateCA(CA);
	SparseMatrix<double> GM(N_plus.size(), N_plus.size());
	SparseMatrix<double> GC_plus(N_plus.size(), N_plus.size());
	SparseMatrix<double> GC_minus(N_plus.size(), N_plus.size());
	GM.setZero();
	GC_plus.setZero();
	GC_minus.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem(mesh,par[1],par[2], GM, GC_plus, GC_minus);
	Solver(GM, GC_plus, GC_minus, Bv,CA);
	cout << "Done reaction!\n";
}