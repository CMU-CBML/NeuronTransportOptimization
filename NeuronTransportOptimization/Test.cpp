#include <stdio.h>
#include <stdlib.h>
#include "Test.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;
const double PI = 3.1415926;
double TEST::calculate_sum(double t, double x,double y, double z)
{
	double sum(0);
	double a, b, c, d;
	double v[3] = { 5,5,5 };
	double diffuseCo=0.01;
	double H(0), Omega(0);
	double n=3;
	a = -(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/ (4 * diffuseCo);
	b = v[0] / (2 * diffuseCo);
	c = v[1] / (2 * diffuseCo);
	d = v[2] / (2 * diffuseCo);

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			for (int k = 1; k <= n; k++)
			{
				H = 8.0 * pow(PI, 3)*i*j*k*(diffuseCo*((1.0 - b)*(1.0 - b) + c*c + d*d) + a)
					*(1.0 + exp(L*(1 - b))*pow(-1.0, i + 1)) / (pow(L*(1.0 - b), 2) + pow(PI*i, 2))
					*(1.0 + exp(-L*c)*pow(-1.0, j + 1)) / (pow(L*c, 2) + pow(PI*j, 2))
					*(1.0 + exp(-L*d)*pow(-1.0, k + 1)) / (pow(L*d, 2) + pow(PI*k, 2));
				Omega = diffuseCo*((PI*i / L)*(PI*i / L) + (PI*j / L)*(PI*j / L) + (PI*k / L)*(PI*k / L));
				sum =sum+ H / (Omega - a)*(exp(-a*t) - exp(-Omega*t))*sin(PI*i*x / L)*sin(PI*j*y / L)*sin(PI*k*z / L);
			}
		}
	}
	//sum *= exp(a*t);
	sum *= exp(a*t + b*x + c*y + d*z);
	sum += exp(x);
	return sum;
}
double TEST::test_function(double x, double y, double z,double &n0,double &n_plus)
{
	double result;
	double A, B, J;
	double lamda = 0.2;
	double q;
	q = sqrt(par[3] / par[0]);
	//A = -(CAi - CAs*exp(-sqrt(par[3] / par[0]))) / (1 - exp(-2.*sqrt(par[3] / par[0])))*sqrt(par[0] / par[3])*(par[3] / par[1]);
	//B = (CAs - CAi*exp(-sqrt(par[3] / par[0]))) / (1 - exp(-2.*sqrt(par[3] / par[0])))*sqrt(par[0] / par[3])*(par[3] / par[1]);
	//J = (CAi - CAs*exp(-sqrt(par[3] / par[0]))) / (1 - exp(-2.*sqrt(par[3] / par[0])))*sqrt(par[0] / par[3])*(par[3] / par[1])*(1 + sqrt(par[3] / par[0])*(par[1] / par[3])) + (lamda - 1)*CAi;
	A = -par[3] * CAi / (q*par[1]);
	B = par[3] * CAs / (q*par[1]);
	J = CAi*(sqrt(par[3] * par[0]) + lamda*par[1]);
	//result = sin(PI*x )*sin(PI*y)*sin(PI*z);
	//result = sin(PI*x / 10.)*sin(PI*y / 10.)*sin(PI*z / 10.);
	//result = sin(PI*x/100.)*sin(PI*y/100.)*sin(PI*z/100.);
	//n0 = A*par[1] / par[3] * (-sqrt(par[3] / par[0]))*exp(-sqrt(par[3] / par[0])*x) + B*par[1] / par[3] * (sqrt(par[3] / par[0]))*exp(sqrt(par[3] / par[0])*(x-L));
	//n_plus = J / par[1] + A*exp(-sqrt(par[3] / par[0])*x) + B*(par[1] / par[3])*sqrt(par[3] / par[0])*exp(sqrt(par[3] / par[0])*(x - L));
	n0 = -par[1] / par[3] * (A*q*exp(-q*x) - B*q*exp(q*(x - L)));
	n_plus = J / par[1] + A*exp(-q*x) + B*exp(q*(x - L));
	return result;
}
TEST::TEST()
{
}

void TEST::InitializeProblem(const vector<double>& CA0, const vector<double>& var, double CAi0, double CAs0, double tstep, double length)
{
	CA.resize(CA0.size());
	CA_Numerical.resize(CA0.size());
	Nplus_Numerical.resize(CA0.size());
	N_plus.resize(CA0.size());
	N_minus.resize(CA0.size());
	error_point.resize(CA0.size());
	Nplus_error.resize(CA0.size());
	par.resize(var.size());//Dn0, v_plus, v_minus// k+, k- //k'+,k'-
	par = var;
	dt = tstep;
	L = length;
	CA = CA0;
	CAi = CAi0;
	CAs = CAs0;
}

void TEST::GaussInfo(int ng)
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
		wght[0] = 1.;			wght[1] = 1.;
		break;
	}
	case 3:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.1127016653792583;			Gpt[1] = 0.5;			Gpt[2] = 0.8872983346207417;
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
		Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;			Gpt[2] = 0.5;			Gpt[3] = 0.7692346550528415;  Gpt[4] = 0.953089922969332;
		wght[0] = 0.2369268850561891;			wght[1] = 0.4786286704993665;			wght[2] = 0.5688888888888889;			wght[3] = 0.4786286704993665; wght[4] = 0.2369268850561891;
		break;
	}
	default:
	{
		Gpt.resize(2);
		wght.resize(2);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;			wght[1] = 1.;
		break;
	}
	}
}

void TEST::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<int>& IEN, double Nx[8], double dNdx[8][3], double& error,double t, double& detJ)
{
	double Nu[2] = { (1. - u), u };
	double Nv[2] = { (1. - v), v };
	double Nw[2] = { (1. - w), w };
	double dNdu[2] = { -1., 1. };
	double dNdv[2] = { -1., 1. };
	double dNdw[2] = { -1., 1. };
	double dNdt[8][3];
	double X(0), Y(0), Z(0);
	double CA_h(0),Nplus_h(0);

	int i, j, k, a, b, loc(0);

	error = 0.;	
	loc = 0;

	for (i = 0; i<2; i++)
	{
		for (k = 0; k < 2; k++)
		{
			Nx[loc] = Nu[k] * Nv[0] * Nw[i];
			X += pt[loc][0] * Nx[loc];
			Y += pt[loc][1] * Nx[loc];
			Z += pt[loc][2] * Nx[loc];
			CA_h += CA_Numerical[IEN[loc]] * Nx[loc];
			Nplus_h += Nplus_Numerical[IEN[loc]] * Nx[loc];
			dNdt[loc][0] = dNdu[k] * Nv[0] * Nw[i];
			dNdt[loc][1] = Nu[k] * dNdv[0] * Nw[i];
			dNdt[loc][2] = Nu[k] * Nv[0] * dNdw[i];
			loc++;
		}
		for (k = 1; k >= 0; k--)
		{
			Nx[loc] = Nu[k] * Nv[1] * Nw[i];
			X += pt[loc][0] * Nx[loc];
			Y += pt[loc][1] * Nx[loc];
			Z += pt[loc][2] * Nx[loc];
			CA_h += CA_Numerical[IEN[loc]] * Nx[loc];
			Nplus_h += Nplus_Numerical[IEN[loc]] * Nx[loc];
			dNdt[loc][0] = dNdu[k] * Nv[1] * Nw[i];
			dNdt[loc][1] = Nu[k] * dNdv[1] * Nw[i];
			dNdt[loc][2] = Nu[k] * Nv[1] * dNdw[i];
			loc++;
		}
	}
	Matrix3d dxdt = Matrix3d::Zero();
	loc = 0;
	for (i = 0; i<2; i++)
	{
		for (j = 0; j<2; j++)
		{
			for (k = 0; k < 2; k++)
			{
				for (a = 0; a<3; a++)
				{
					for (b = 0; b<3; b++)
					{
						dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
					}
				}
				loc++;
			}
		}
	}

	Matrix3d dtdx = dxdt.inverse();
	detJ = dxdt.determinant();
	detJ = 0.125*detJ;
	error = pow((CA_h - calculate_sum(t, X, Y, Z)), 2);
//	error = (CA_h - test_function(X, Y, Z))*(CA_h - test_function(X, Y, Z));

}

void TEST::errorL2(const vector<Element3D>& mesh, double &result,double t)
{
	double ERROR_L2(0.);
	double er,element_error;
	unsigned int e, i, j, k;
	double Nx[8];
	double dNdx[8][3];
	double detJ;
	GaussInfo(2);
	
	for (e = 0; e < mesh.size(); e++)
	{
		element_error = 0.;
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], mesh[e].pts, mesh[e].IEN, Nx, dNdx, er,t, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					element_error += er*detJ;
				}
			}
		}
		ERROR_L2 += element_error;
	}
	cout <<"Error Calculation Done!\n";
	result=sqrt(ERROR_L2);
}
void TEST::Run(const vector<array<double, 3>>& pts, const vector<Element3D>& bzmesh, const vector<double>& Bv, vector<double>& CA1, const vector<int>& pid_loc, const double time, const int judge, string fn)
{

	if (judge == 1)//cube
	{
		for (int i = 0; i < pts.size(); i++)
		{
			//CA[i] = CAs + (CAi - CAs)*calculate_sum(time,pts[i][0], pts[i][1], pts[i][2]); //cube
			//test_function(pts[i][0], pts[i][1], pts[i][2],CA[i],N_plus[i]);
			CA[i] = calculate_sum(time, pts[i][0], pts[i][1], pts[i][2]);
		}
	}

}
void TEST::VisualizeVTK(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn)
{
	string fname = fn + "_N0_analytical.vtk";
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
		fout << "POINT_DATA " << CA.size() << "\nSCALARS N0 float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<CA.size(); i++)
		{
			fout << CA[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//fname = fn + "_Nplus_analytical.vtk";
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << spt.size() << " float\n";
	//	for (i = 0; i<spt.size(); i++)
	//	{
	//		fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
	//	}
	//	fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
	//	for (i = 0; i<mesh.size(); i++)
	//	{
	//		fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
	//			<< " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
	//	}
	//	fout << "\nCELL_TYPES " << mesh.size() << '\n';
	//	for (i = 0; i<mesh.size(); i++)
	//	{
	//		fout << "12\n";
	//	}
	//	fout << "POINT_DATA " << N_plus.size() << "\nSCALARS N+ float 1\nLOOKUP_TABLE default\n";
	//	for (uint i = 0; i<N_plus.size(); i++)
	//	{
	//		fout << N_plus[i] << "\n";
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}
}

void TEST::VisualizeVTK_ERROR(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn)
{
	string fname = fn + "_N0_error.vtk";
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
		fout << "POINT_DATA " << CA.size() << "\nSCALARS Nplus_error float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<CA.size(); i++)
		{
			fout << error_point[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
	//fname = fn + "_Nplus_error.vtk";
	//fout.open(fname.c_str());
	//if (fout.is_open())
	//{
	//	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout << "POINTS " << spt.size() << " float\n";
	//	for (i = 0; i<spt.size(); i++)
	//	{
	//		fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
	//	}
	//	fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
	//	for (i = 0; i<mesh.size(); i++)
	//	{
	//		fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
	//			<< " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
	//	}
	//	fout << "\nCELL_TYPES " << mesh.size() << '\n';
	//	for (i = 0; i<mesh.size(); i++)
	//	{
	//		fout << "12\n";
	//	}
	//	fout << "POINT_DATA " << Nplus_error.size() << "\nSCALARS N0_error float 1\nLOOKUP_TABLE default\n";
	//	for (uint i = 0; i<CA.size(); i++)
	//	{
	//		fout << Nplus_error[i] << "\n";
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}
}