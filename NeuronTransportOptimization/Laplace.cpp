#include "Laplace.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;
const double PI = 3.1415926;



//Poisson Test
double Laplace::test_function(double x, double y, double z)
{
	double result;
	result = sin(PI*x/100.)*sin(PI*y/100.)*sin(PI*z/100.);
	return result;
}

Laplace::Laplace()
{
}

vector<double> Laplace::GetCA()
{
	return CA;
}

vector<double> Laplace::GetNplus()
{
	return N_plus;
}

vector<double> Laplace::GetNminus()
{
	return N_minus;
}

void Laplace::GaussInfo(int ng)
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
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;							wght[1]=1.;
			break;
		}
	}
}

void Laplace::InitializeProblem(vector<Element3D> &bzmesh, vector<array<double, 3>> &velocity_node, const vector<double>& CA0, const vector<double>& N_plus0, const vector<double>& N_minus0, const vector<double>& var, double tstep)
{
	int i, j, e;
	
	cout << "Initializing Problem..." << endl;

	shape.Nx.resize(bzpt_num);
	shape.dNdx.resize(bzpt_num);

	CA.resize(CA0.size());
	N_plus.resize(N_plus0.size());
	N_minus.resize(N_minus0.size());
	vplus.resize(CA0.size());
	vminus.resize(CA0.size());
	MatrixSolve.resize(CA0.size() * 2, CA0.size() * 2);
	VectorSolve.resize(CA0.size() * 2, 0);

	if (var.size() != 0)
	{
		DA = var[0];
		par = var;
	}
	else
	{
		cerr << "0 variables!\n"; getchar();
	}
	dt = tstep;

	judge = 0;
	CA = CA0;
	N_plus = N_plus0;
	N_minus = N_minus0;

	//for (e = 0; e < bzmesh.size(); e++)
	//{
	//	Tmesh2Bzmesh_velocity(bzmesh[e], vplus, bzmesh[e].vplus_bz);
	//	Tmesh2Bzmesh_velocity(bzmesh[e], vminus, bzmesh[e].vminus_bz);
	//	Tmesh2Bzmesh_value(bzmesh[e], CA, bzmesh[e].CA_bz);
	//	Tmesh2Bzmesh_value(bzmesh[e], N_plus, bzmesh[e].Nplus_bz);
	//	Tmesh2Bzmesh_value(bzmesh[e], N_minus, bzmesh[e].Nminus_bz);
	//}
}

void Laplace::CoordinateTransfer(double u, double v, double w, double &x, double &y, double &z, const vector<array<double, 3>>& pt)
{
	int i, j, k, l, loc(0);
	x = 0;	y = 0; z = 0;
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				x += pt[loc][0] * Nu[k] * Nv[j] * Nw[i];
				y += pt[loc][1] * Nu[k] * Nv[j] * Nw[i];
				z += pt[loc][2] * Nu[k] * Nv[j] * Nw[i];
				loc++;
			}
		}
	}	
}

void Laplace::VelocityTransfer(double u, double v, double w, const vector<array<double,64>> &cmat, const vector<array<double, 3>>& v_node, double v_tmp[3])
{
	int i, j, k, l, loc(0);
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	for (l = 0; l < dim; l++)
	{
		v_tmp[l] = 0;
	}

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (l = 0; l < dim; l++)
				{
					v_tmp[l] += v_node[loc][l] * Nu[k] * Nv[j] * Nw[i];
					loc++;
				}
			}
		}
	}
}

void Laplace::Bzmesh2Tmesh_value(Element3D &bzel, const double value_bzpt[bzpt_num], vector<double>& value_cpt)
{
	int i, j, k;

	for (i = 0; i < bzpt_num; i++)
	{
		for (j = 0; j < bzel.IEN.size(); j++)
		{
			value_cpt[bzel.IEN[j]] += bzel.cmat[j][i] * value_bzpt[i];
		}
	}
}

void Laplace::Tmesh2Bzmesh_value(Element3D &bzel, const vector<double>& value_cpt, double value_bzpt[bzpt_num])
{
	int i, j, k;
	for (i = 0; i < bzpt_num; i++)
			value_bzpt[i] = 0.;

	for (i = 0; i < bzpt_num; i++)
	{
		for (j = 0; j < bzel.IEN.size(); j++)
		{
			value_bzpt[i] += bzel.cmat[j][i] * value_cpt[bzel.IEN[j]];
		}
	}
}

void Laplace::Tmesh2Bzmesh_velocity(Element3D &bzel, const vector<array<double, 3>>& v_cpt_node, vector<array<double, 3>>& v_bzpt_node)
{
	int i, j, k;
	for (i = 0; i < bzpt_num; i++)
		for(j=0;j<3;j++)
			v_bzpt_node[i][j] = 0.;

	for (i = 0; i < bzpt_num; i++)
	{
		for (j = 0; j < bzel.IEN.size(); j++)
		{
			for (k = 0; k < 3; k++)
			{
				v_bzpt_node[i][k] += bzel.cmat[j][i] * v_cpt_node[bzel.IEN[j]][k];
			}
		}		
	}
}

void Laplace::SUPGcoefficient(double velocity[3], double s, const double dudx[3][3], double &tau)
{
	double b[3] = {0, 0, 0};
	double Pe, Da, ksi, h;
	double c1(1.), c2(2.5);
	double norm_v = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b[i] += velocity[j] * dudx[i][j];
		}
	}

	double norm_b = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
	
	h = 2 * norm_v / norm_b;
	Pe = 0.5*norm_v*h / par[0];
	ksi = c2*(1 / tanh(c1*Pe / c2) - (c2 / (c1*Pe)));

	tau = 0.5*norm_v*h*ksi;
	tau = 0.5*norm_v*h*c2;
	tau = tau / (norm_v*norm_v);

}

void Laplace::SUPGcoefficient2(double v[3], double dNdx[bzpt_num][dim], double &tau_supg)
{
	double tau1(0.), tau2(0.), tau3(0.);

	for (int i = 0; i < bzpt_num; i++)
	{
		tau1 += abs(v[0] * dNdx[i][0] + v[1] * dNdx[i][1] + v[2] * dNdx[i][2]);
	}
	tau1 = 1. / tau1;
	tau2 = dt / 2.;
	tau3 = 0.;

	//tau_supg = 1. / sqrt(1. / (tau1*tau1) + 1. / (tau2*tau2) + 1. / (tau3*tau3));
	tau_supg = 1. / sqrt(1. / (tau1*tau1) + 1. / (tau2*tau2));
}

void Laplace::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[bzpt_num], double dNdx[bzpt_num][dim], double dudx[3][3], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[64][3];
	int i, j, k, a, b, loc(0);
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nx[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	Matrix3d dxdt = Matrix3d::Zero();
	for(loc=0;loc<bzpt_num;loc++)
	{
		for (a = 0; a<3; a++)
		{
			for (b = 0; b<3; b++)
			{
				dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
			}
		}

	}
	//cout << dxdt << "\n\n";
	Matrix3d dtdx = dxdt.inverse();
	for (i = 0; i<bzpt_num; i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			dudx[i][j] = dtdx(i, j);

	detJ = dxdt.determinant();
	detJ = 0.125*detJ;
}

void Laplace::WeightingFunction(const double velocity[3], const double& tau, const double Nx[bzpt_num], const double dNdx[bzpt_num][dim], double Wx[bzpt_num])
{
	double U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
	for (int loc = 0; loc < bzpt_num; loc++)
	{
		//Nx[loc] += k_bar / (U*U)*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
		//Nx[loc] += k_bar[loc] / (U*U)*(velocity[loc][0] * dNdx[loc][0] + velocity[loc][1] * dNdx[loc][1] + velocity[loc][2] * dNdx[loc][2]);
		Wx[loc] = Nx[loc] + tau*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
	}
}

void Laplace::ElementVelocity(double Nx[bzpt_num], const vector<array<double, 3>> v_node, double v_tmp[3])
{
	for (int j = 0; j < dim; j++)
	{
		v_tmp[j] = 0;
	}
	for (int i = 0; i < bzpt_num; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			v_tmp[j] += v_node[i][j] * Nx[i];
		}
	}
}

void Laplace::ElementValue(double Nx[bzpt_num], const double value_bz[bzpt_num], double &value)
{
	value = 0.;
	for (int i = 0; i < bzpt_num; i++)
	{
		value += value_bz[i] * Nx[i];
	}
}

void Laplace::ElementMassMatrix(double Weight[bzpt_num], double Nx[bzpt_num], double detJ, double EM[bzpt_num][bzpt_num])
{
	int i, j;
	for (i = 0; i<bzpt_num; i++)
	{
		for (j = 0; j<bzpt_num; j++)
		{
			EM[i][j] = Weight[i] * Nx[j] * detJ;
		}
	}
}

void Laplace::ElementConvectionMatrix(double Nx[bzpt_num], double dNdx[bzpt_num][dim], double v[dim], double detJ, double EC[bzpt_num][bzpt_num])
{
	int i, j;
	for (i = 0; i<bzpt_num; i++)
	{
		for (j = 0; j<bzpt_num; j++)
		{
			EC[i][j] = Nx[i] * (v[0] * dNdx[j][0] + v[1] * dNdx[j][1] + v[2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementStiffMatrix(double dNdx[bzpt_num][dim], double detJ, double EK[bzpt_num][bzpt_num])
{
	int i, j;
	for (i = 0; i<bzpt_num; i++)
	{
		for (j = 0; j<bzpt_num; j++)
		{
			EK[i][j] = (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementLoadVector(double Fx[bzpt_num], double detJ, double EF[bzpt_num])
{
	int i;
	for (i = 0; i<bzpt_num; i++)
	{
		EF[i] = Fx[i] * detJ;
	}
}

void Laplace::Tangent(double EM[bzpt_num][bzpt_num], double EM_plus[bzpt_num][bzpt_num], double EM_minus[bzpt_num][bzpt_num], double EK[bzpt_num][bzpt_num], double EC_plus[bzpt_num][bzpt_num], double EC_minus[bzpt_num][bzpt_num], double EMatrixSolve[bzpt_num * 2][bzpt_num * 2])
{
	int i, j, A, B;
	for (i = 0; i < bzpt_num; i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			EMatrixSolve[i + 0][j + 0] += EM[i][j] + dt*DA*EK[i][j] + dt*(par[3] + par[4]) * EM[i][j];
			EMatrixSolve[i + 0][j + bzpt_num] += -par[5] * dt * EM[i][j];
			EMatrixSolve[i + bzpt_num][j + 0] += -par[3] * dt * EM_plus[i][j];
			EMatrixSolve[i + bzpt_num][j + bzpt_num] += EM_plus[i][j] + dt*EC_plus[i][j] + dt*par[5] * EM_plus[i][j];

			//EMatrixSolve[i + 0][j + 16] += -par[6] * dt * EM[i][j];
			//EMatrixSolve[i + 16][j + 0] += -par[4] * dt * EM_minus[i][j];
			//EMatrixSolve[i + 16][j + 16] += EM_minus[i][j] + dt*EC_minus[i][j] + dt*par[6] * EM_minus[i][j];
		}
	}
	
}

void Laplace::Residual(const double CA_bz[bzpt_num], const double Nplus_bz[bzpt_num], const double Nminus_bz[bzpt_num], double EM[bzpt_num][bzpt_num], double EM_plus[bzpt_num][bzpt_num], double EM_minus[bzpt_num][bzpt_num], double EVectorSolve[bzpt_num * 2])
{
	int i, j, A;
	for (i = 0; i < bzpt_num; i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			EVectorSolve[i] += EM[i][j] * CA_bz[j];
			EVectorSolve[i + bzpt_num] += EM_plus[i][j] * Nplus_bz[j];
			//EVectorSolve[i + bzpt_num*2] += EM_minus[i][j] * Nminus_bz[j];
		}
	}

}

void Laplace::BuildLinearSystem(const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<array<double, 3>> &cpts, const vector<int>& pid_loc, const vector<int>& label, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus, SparseMatrix<double, RowMajor>& MatrixSolve, vector<double>& VectorSolve)
{
	unsigned int e, i, j, k, m;
	double detJ;
	double tau[2];
	double dudx[3][3];

	double EM[bzpt_num][bzpt_num], EK[bzpt_num][bzpt_num], EC_plus[bzpt_num][bzpt_num], EC_minus[bzpt_num][bzpt_num], EM_plus[bzpt_num][bzpt_num], EM_minus[bzpt_num][bzpt_num];
	double EMatrixSolve[bzpt_num * 2][bzpt_num * 2];
	double Fx[bzpt_num], EF[bzpt_num], EVectorSolve[bzpt_num * 2];
	double Nx[bzpt_num], Npx[bzpt_num], Nmx[bzpt_num];
	double dNdx[bzpt_num][dim];
	vector<vector<double>> EMatrixSolve1;
	vector<double>EVectorSolve1;

	//vector<vector<double>> EM, EK, EC_plus, EC_minus, EM_plus, EM_minus;
	//vector<vector<double>> EMatrixSolve;
	//vector<double>Fx, EF, EVectorSolve;
	//vector<double>Nx, Npx, Nmx;
	//vector<array<double, 3>> dNdx;

	vector<Triplet<double>> trilist_MatrixSolve;
	for (e = 0; e<bzmesh.size(); e++)
	{
		//cout << "Element: " << e << endl;
		int nen = bzmesh[e].IEN.size();
		EMatrixSolve1.clear(); EMatrixSolve1.resize(nen * 2);
		for (i = 0; i < nen * 2; i++)
			EMatrixSolve1[i].resize(nen * 2, 0.0);
		EVectorSolve1.clear(); EVectorSolve1.resize(nen * 2, 0.0);

		//EM.clear();		EK.clear();		EC_plus.clear();	EC_minus.clear();	EM_plus.clear();	EM_minus.clear();
		//EMatrixSolve.clear();
		//Fx.clear(); EF.clear(); EVectorSolve.clear();
		//Nx.clear(); Npx.clear(); Nmx.clear();
		//dNdx.clear();
		//
		//EM.resize(nen);		EK.resize(nen);		EC_plus.resize(nen);	EC_minus.resize(nen);	EM_plus.resize(nen);	EM_minus.resize(nen);
		//EMatrixSolve.resize(nen * 2);
		//Fx.resize(nen); EF.resize(nen); EVectorSolve.resize(nen*2);
		//Nx.resize(nen); Npx.resize(nen); Nmx.resize(nen); 
		//dNdx.resize(nen);
		//
		//for (i = 0; i < nen; i++)
		//{
		//	EM[i].resize(nen);		EK[i].resize(nen);		EC_plus[i].resize(nen);		EC_minus[i].resize(nen);	EM_plus[i].resize(nen);	EM_minus[i].resize(nen);
		//	EMatrixSolve[i].resize(nen * 2);
		//	EMatrixSolve[i+nen].resize(nen * 2);
		//}

		double v_plus[3], v_minus[3];
		double v_tmp_plus[3], v_tmp_minus[3];
		double CA_tmp, Nplus_tmp, Nminus_tmp;

		//vector<array<double, 3>> v_cpt_plus, v_cpt_minus;
		//vector<array<double, 3>> v_bzpt_plus, v_bzpt_minus;
		//v_cpt_plus.resize(nen); v_cpt_minus.resize(nen);
		//v_bzpt_plus.resize(bzpt_num); v_bzpt_minus.resize(bzpt_num);

		//vector<double> CA_cpt, Nplus_cpt, Nminus_cpt;
		//double CA_bzpt[bzpt_num], Nplus_bzpt[bzpt_num], Nminus_bzpt[bzpt_num];
		//CA_cpt.resize(nen); Nplus_cpt.resize(nen); Nminus_cpt.resize(nen);

		//for (i = 0; i<bzpt_num; i++)
		//{
		//	for (j = 0; j<bzpt_num; j++)
		//	{
		//		EM[i][j] = 0.;
		//		EK[i][j] = 0.;
		//		EC_plus[i][j] = 0.;
		//		EC_minus[i][j] = 0.;
		//		EM_plus[i][j] = 0.;
		//		EM_minus[i][j] = 0.;
		//	}
		//	EF[i] = 0.;
		//}

		for (i = 0; i < bzpt_num * 2; i++)
		{
			for (j = 0; j < bzpt_num * 2; j++)
			{
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}


		//for (i = 0; i < nen; i++)
		//{
		//	CA_cpt[i] = CA[bzmesh[e].IEN[i]];
		//	Nplus_cpt[i] = N_plus[bzmesh[e].IEN[i]];
		//	Nminus_cpt[i] = N_minus[bzmesh[e].IEN[i]];
		//	for (j = 0; j < dim; j++)
		//	{
		//		v_cpt_plus[i][j] = Vplus*velocity_node[bzmesh[e].IEN[i]][j];
		//		v_cpt_minus[i][j] = Vminus*velocity_node[bzmesh[e].IEN[i]][j];
		//	}
		//}
		//
		//Tmesh2Bzmesh_velocity(bzmesh[e], v_cpt_plus, v_bzpt_plus);
		//Tmesh2Bzmesh_velocity(bzmesh[e], v_cpt_minus, v_bzpt_minus);
		//Tmesh2Bzmesh_value(bzmesh[e], CA_cpt, CA_bzpt);
		//Tmesh2Bzmesh_value(bzmesh[e], Nplus_cpt, Nplus_bzpt);
		//Tmesh2Bzmesh_value(bzmesh[e], Nminus_cpt, Nminus_bzpt);

		//cout << "Gauss Quadrature " << e << endl;

		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, dNdx, dudx, detJ);
					ElementVelocity(Nx, bzmesh[e].vplus_bz, v_tmp_plus);
					ElementVelocity(Nx, bzmesh[e].vminus_bz, v_tmp_minus);
					//ElementValue(Nx, bzmesh[e].CA_bz, CA_tmp);
					//ElementValue(Nx, bzmesh[e].Nplus_bz, Nplus_tmp);
					//ElementValue(Nx, bzmesh[e].Nminus_bz, Nminus_tmp);
					//SUPGcoefficient(v_tmp_plus, par[5], dudx, k_bar[0]);
					//SUPGcoefficient(v_tmp_minus, par[6], dudx, k_bar[1]);
					SUPGcoefficient2(v_tmp_plus, dNdx, tau[0]);
					SUPGcoefficient2(v_tmp_minus, dNdx, tau[1]);

					WeightingFunction(v_tmp_plus, tau[0], Nx, dNdx, Npx);
					WeightingFunction(v_tmp_minus, tau[1], Nx, dNdx, Nmx);
					detJ = wght[i] * wght[j] * wght[k] * detJ;


					ElementMassMatrix(Nx, Nx, detJ, EM);
					ElementMassMatrix(Npx, Nx, detJ, EM_plus);
					ElementMassMatrix(Nmx, Nx, detJ, EM_minus);
					ElementStiffMatrix(dNdx, detJ, EK);
					ElementConvectionMatrix(Npx, dNdx, v_tmp_plus, detJ, EC_plus);
					ElementConvectionMatrix(Nmx, dNdx, v_tmp_minus, detJ, EC_minus);
					Tangent(EM, EM_plus, EM_minus, EK, EC_plus, EC_minus, EMatrixSolve);
					Residual(bzmesh[e].CA_bz, bzmesh[e].Nplus_bz, bzmesh[e].Nminus_bz, EM, EM_plus, EM_minus, EVectorSolve);
				}
			}
		}
		Bzmesh2Tmesh_Matrix(bzmesh[e].cmat, bzmesh[e].IEN, EMatrixSolve, EMatrixSolve1);
		Bzmesh2Tmesh_Vector(bzmesh[e].cmat, bzmesh[e].IEN, EVectorSolve, EVectorSolve1);

		//Apply BCs
		double CA_bc, Nplus_bc, Nminus_bc;
		for (i = 0; i < nen; i++)
		{
			int A = bzmesh[e].IEN[i];
			//Cylinder BCs			
			if (sqrt(pow(cpts[A][1], 2) + pow(cpts[A][2], 2)) > 0.48)
			{
				CA_bc = 0;
				ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
				Nplus_bc = 0;
				ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
				//Nminus_bc = 0;
				//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			}
			else if (label[A] == 1)
			{
				CA_bc = 1.0;
				ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
				Nplus_bc = 2.0;
				ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			}

			////Bifurcation BCs
			//if (pid_loc[A] != -1)
			//{
			//	double R0 = sqrt(pow(cpts[A][0] - 0.0, 2) + pow(cpts[A][1] - 0.0, 2) + pow(cpts[A][2] - 0.0, 2));
			//	double R1 = sqrt(pow(cpts[A][0] - 7.33177, 2) + pow(cpts[A][1] - (-6.29269), 2) + pow(cpts[A][2] - (-5.655), 2));
			//	double R2 = sqrt(pow(cpts[A][0] - 4.85463, 2) + pow(cpts[A][1] - (-5.84302), 2) + pow(cpts[A][2] - 1.701, 2));
			//
			//	if (label[A] == 1)
			//	{
			//		if (R0 > 0.7)//inlet wall
			//		{
			//			CA_bc = 0.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//			Nplus_bc = 0.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//		else//inlet
			//		{
			//			CA_bc = 1.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//			Nplus_bc = 2.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//
			//	}
			//	if (label[A] < 1) //wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 2 && R1>0.45) //outlet1 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 3 && R2>0.47) //outlet2 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//
			//}

			//3 Bifurcation BCs
			//if (pid_loc[A] != -1)
			//{
			//	double R0 = sqrt(pow(cpts[A][0] - 0.0, 2)     + pow(cpts[A][1] - 0.0, 2)     + pow(cpts[A][2] - 0.0, 2));
			//	double R1 = sqrt(pow(cpts[A][0] - 63.6235, 2) + pow(cpts[A][1] - 14.7335, 2) + pow(cpts[A][2] - 8.25363, 2));
			//	double R2 = sqrt(pow(cpts[A][0] - 46.7595, 2) + pow(cpts[A][1] - 22.4332, 2) + pow(cpts[A][2] - 7.46334, 2));
			//	double R3 = sqrt(pow(cpts[A][0] - 35.7126, 2) + pow(cpts[A][1] - 36.7129, 2) + pow(cpts[A][2] - 2.66824, 2));
			//	double R4 = sqrt(pow(cpts[A][0] - 28.7469, 2) + pow(cpts[A][1] - 43.6465, 2) + pow(cpts[A][2] - 2.81018, 2));
			//	
			//	if (label[A] == 1)
			//	{
			//		if (R0 > 2.0)//inlet wall
			//		{
			//			CA_bc = 0.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//			Nplus_bc = 0.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//		else//inlet
			//		{
			//			CA_bc = 1.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//			Nplus_bc = 2.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//
			//	}
			//	if (label[A] < 1) //wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 2 && R1>0.78) //outlet1 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 3 && R2>0.61) //outlet2 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 4 && R3>0.64) //outlet3 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 5 && R4>0.71) //outlet4 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve1, EVectorSolve1);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve1, EVectorSolve1);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//
			//}
		}
		Assembly(EMatrixSolve1, EVectorSolve1, bzmesh[e].IEN, trilist_MatrixSolve, VectorSolve);
	}
	//cout << "Assemble Global Matrix" << endl;
	if (judge == 0)
		MatrixSolve.setFromTriplets(trilist_MatrixSolve.begin(), trilist_MatrixSolve.end());
}


//Using Control Points Data Structure
void Laplace::SUPGcoefficientIGA(double s, double v[3], vector<array<double, 3>>& dNdx, double &tau_supg)
{
	double tau1(0.), tau2(0.), tau3(0.);

	for (int i = 0; i < dNdx.size(); i++)
	{
		tau1 += abs(v[0] * dNdx[i][0] + v[1] * dNdx[i][1] + v[2] * dNdx[i][2]);
	}
	tau1 = 1. / tau1;
	tau2 = dt / 2.;
	tau3 = 0.;

	//tau_supg = 1. / sqrt(1. / (tau1*tau1) + 1. / (tau2*tau2) + 1. / (tau3*tau3));
	tau_supg = 1. / sqrt(1. / (tau1*tau1) + 1. / (tau2*tau2) + 1. / (s*s));
}

void Laplace::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, double dudx[3][3], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[bzpt_num][3];
	double Nx_bz[bzpt_num];
	double dNdx_bz[bzpt_num][3];

	Nx.clear();
	dNdx.clear();
	Nx.resize(cmat.size(), 0);
	dNdx.resize(cmat.size(), { 0, 0, 0 });

	int i, j, k, a, b, c, loc;
	loc = 0;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}

	Matrix3d dxdt = Matrix3d::Zero();
	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (a = 0; a<3; a++)
		{
			for (b = 0; b<3; b++)
			{
				dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
			}
		}
	}
	Matrix3d dtdx = dxdt.inverse();

	//1st derivatives
	for (i = 0; i<bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx_bz[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx_bz[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			dudx[i][j] = dtdx(i, j);
		}
	}
	detJ = dxdt.determinant();
	detJ = 0.125*detJ;

	for (i = 0; i < cmat.size(); i++)
	{
		for (j = 0; j < bzpt_num; j++)
		{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 3; m++)
			{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
			}
		}
	}
}

void Laplace::WeightingFunction(const double velocity[3], const double& tau, const vector<double> &Nx, const vector<array<double, 3>> &dNdx, vector<double> &Wx)
{
	double U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
	for (int loc = 0; loc < Nx.size(); loc++)
	{
		//Nx[loc] += k_bar / (U*U)*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
		//Nx[loc] += k_bar[loc] / (U*U)*(velocity[loc][0] * dNdx[loc][0] + velocity[loc][1] * dNdx[loc][1] + velocity[loc][2] * dNdx[loc][2]);
		Wx[loc] = Nx[loc] + tau*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
	}
}

void Laplace::ElementValue(const vector<double> &Nx, const vector<double> value_node, double &value)
{
	value = 0.;
	for (int i = 0; i < Nx.size(); i++)
	{
		value += value_node[i] * Nx[i];
	}
}

void Laplace::ElementVelocity(const vector<double> &Nx, const vector<array<double, 3>> &v_node, double v_tmp[3])
{
	for (int j = 0; j < dim; j++)
	{
		v_tmp[j] = 0;
	}
	for (int i = 0; i < Nx.size(); i++)
	{
		for (int j = 0; j < dim; j++)
		{
			v_tmp[j] += v_node[i][j] * Nx[i];
		}
	}
}

void Laplace::ElementMassMatrix(vector<double> &Weight, vector<double>& Nx, double detJ, vector<vector<double>>& EM)
{
	int i, j;
	for (i = 0; i<Nx.size(); i++)
	{
		for (j = 0; j<Nx.size(); j++)
		{
			EM[i][j] = Weight[i] * Nx[j] * detJ;
		}
	}
}

void Laplace::ElementConvectionMatrix(vector<double>& Nx, vector<array<double, 3>>& dNdx, double v[3], double detJ, vector<vector<double>>& EC)
{
	int i, j;
	for (i = 0; i<Nx.size(); i++)
	{
		for (j = 0; j<Nx.size(); j++)
		{
			EC[i][j] = Nx[i] * (v[0] * dNdx[j][0] + v[1] * dNdx[j][1] + v[2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementStiffMatrix(vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK)
{
	int i, j;
	for (i = 0; i<dNdx.size(); i++)
	{
		for (j = 0; j<dNdx.size(); j++)
		{
			EK[i][j] = (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::Tangent(const int nen, vector<vector<double>>& EM, vector<vector<double>>& EM_plus, vector<vector<double>>& EM_minus, vector<vector<double>>& EK, vector<vector<double>>& EC_plus, vector<vector<double>>& EC_minus, vector<vector<double>>& EMatrixSolve)
{
	int i, j, A, B;
	for (i = 0; i < nen; i++)
	{
		for (j = 0; j < nen; j++)
		{
			EMatrixSolve[i + 0][j + 0] += EM[i][j] + dt*DA*EK[i][j] + dt*(par[3] + par[4]) * EM[i][j];
			EMatrixSolve[i + 0][j + nen] += -par[5] * dt * EM[i][j];
			EMatrixSolve[i + nen][j + 0] += -par[3] * dt * EM_plus[i][j];
			EMatrixSolve[i + nen][j + nen] += EM_plus[i][j] + dt*EC_plus[i][j] + dt*par[5] * EM_plus[i][j];

			//EMatrixSolve[i + 0][j + 16] += -par[6] * dt * EM[i][j];
			//EMatrixSolve[i + 16][j + 0] += -par[4] * dt * EM_minus[i][j];
			//EMatrixSolve[i + 16][j + 16] += EM_minus[i][j] + dt*EC_minus[i][j] + dt*par[6] * EM_minus[i][j];
		}
	}
}

void Laplace::Residual(const int nen, const double CA, const double Nplus, const double Nminus, const vector<double> &Nx, const vector<double> &Npx, const vector<double> &Nmx, const double detJ, vector<double> &EVectorSolve)
{
	int i, j, A;
	for (i = 0; i < nen; i++)
	{
		EVectorSolve[i] += CA* Nx[i] * detJ;
		EVectorSolve[i + nen] += Nplus* Npx[i] * detJ;
		//EVectorSolve[i + nen * 2] += Nminus*Nmx[i] * detJ;
	}
}

void Laplace::Bzmesh2Tmesh_Matrix(const vector<array<double, 64>>& cmat, const vector<int>& IEN, double EMatrixSolve[bzpt_num * 2][bzpt_num * 2], vector<vector<double>> &EMatrixSolve1)
{
	int i, j, A, B;
	int nen = IEN.size();
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			for (A = 0; A<bzpt_num; A++)
			{
				for (B = 0; B<bzpt_num; B++)
				{
					EMatrixSolve1[i + 0][j + 0] += cmat[i][A] * cmat[j][B] * EMatrixSolve[A + 0][B + 0];
					EMatrixSolve1[i + 0][j + nen] += cmat[i][A] * cmat[j][B] * EMatrixSolve[A + 0][B + bzpt_num];
					EMatrixSolve1[i + nen][j + 0] += cmat[i][A] * cmat[j][B] * EMatrixSolve[A + bzpt_num][B + 0];
					EMatrixSolve1[i + nen][j + nen] += cmat[i][A] * cmat[j][B] * EMatrixSolve[A + bzpt_num][B + bzpt_num];
				}
			}
		}
	}
}

void Laplace::Bzmesh2Tmesh_Vector(const vector<array<double, 64>>& cmat, const vector<int>& IEN, double EVectorSolve[bzpt_num * 2], vector<double> &EVectorSolve1)
{
	int i, A;
	int nen = IEN.size();
	for (i = 0; i < IEN.size(); i++)
	{
		for (A = 0; A < bzpt_num; A++)
		{
			EVectorSolve1[i] += cmat[i][A] * EVectorSolve[A];
			EVectorSolve1[i + nen] += cmat[i][A] * EVectorSolve[A + bzpt_num];
		}
	}
}

//void Laplace::BasisFunction(double u, double v, double w, const int nen, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, vector<double> &Fx,double dudx[3][3], double& detJ)
//{
//	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
//	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
//	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
//	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
//	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
//	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
//	double dNdt[bzpt_num][3];
//	double Nx_bz[bzpt_num];
//	double dNdx_bz[bzpt_num][3];
//	double X(0), Y(0), Z(0);
//
//	int i, j, k, a, b, loc(0);
//	Nx.clear();
//	dNdx.clear();
//	Nx.resize(nen,0);
//	dNdx.resize(nen, {0});
//	loc = 0;
//	for (i = 0; i<4; i++)
//	{
//		for (j = 0; j<4; j++)
//		{
//			for (k = 0; k < 4; k++)
//			{
//				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
//				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
//				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
//				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
//				loc++;
//			}
//		}
//	}
//	Matrix3d dxdt = Matrix3d::Zero();
//	loc = 0;
//	for(loc=0;loc<bzpt_num;loc++)
//	{
//		for (a = 0; a<3; a++)
//		{
//			for (b = 0; b<3; b++)
//			{
//				dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
//			}
//		}
//	}
//	
//	Matrix3d dtdx = dxdt.inverse();
//	for (i = 0; i<bzpt_num; i++)
//	{
//		dNdx_bz[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
//		dNdx_bz[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
//		dNdx_bz[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
//	}
//
//	for (i = 0; i < 3; i++)
//		for (j = 0; j < 3; j++)
//			dudx[i][j] = dtdx(i, j);
//
//	detJ = dxdt.determinant();
//	detJ = 0.125*detJ;
//
//	for (i = 0; i < nen; i++)
//	{
//		for (j = 0; j < bzpt_num; j++)
//		{
//			Nx[i] += cmat[i][j] * Nx_bz[j];
//			for (int m = 0; m < 3; m++)
//			{
//				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
//			}
//		}
//	}
//}
//
//void Laplace::WeightingFunction(double u, double v, double w, const int nen, const double velocity[3],const double& k_bar, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, double& detJ)
//{
//	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
//	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
//	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
//	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
//	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
//	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
//	double dNdt[bzpt_num][3];
//	double X(0), Y(0), Z(0);
//	
//
//	double Nx_bz[bzpt_num];
//	double dNdx_bz[bzpt_num][3];
//
//	int i, j, k, a, b, loc(0);
//	Nx.clear();
//	dNdx.clear();
//	Nx.resize(nen, 0);
//	dNdx.resize(nen, { 0 });
//
//	loc = 0;
//	for (i = 0; i<4; i++)
//	{
//		for (j = 0; j<4; j++)
//		{
//			for (k = 0; k < 4; k++)
//			{
//				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
//				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
//				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
//				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
//				loc++;
//			}
//		}
//	}
//
//	Matrix3d dxdt = Matrix3d::Zero();
//	loc = 0;
//	for (i = 0; i<4; i++)
//	{
//		for (j = 0; j<4; j++)
//		{
//			for (k = 0; k < 4; k++)
//			{
//				for (a = 0; a<3; a++)
//				{
//					for (b = 0; b<3; b++)
//					{
//						dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
//					}
//				}
//				loc++;
//			}
//		}
//	}
//
//	Matrix3d dtdx = dxdt.inverse();
//	for (i = 0; i<bzpt_num; i++)
//	{
//		dNdx_bz[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
//		dNdx_bz[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
//		dNdx_bz[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
//	}
//
//	for (i = 0; i < nen; i++)
//	{
//		for (j = 0; j < bzpt_num; j++)
//		{
//			Nx[i] += cmat[i][j] * Nx_bz[j];
//			for (int m = 0; m < 3; m++)
//			{
//				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
//			}
//		}
//	}
//
//	double U = sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
//	for (loc = 0; loc < nen; loc++)
//	{
//		//Nx[loc] += k_bar / (U*U)*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
//		//Nx[loc] += k_bar[loc] / (U*U)*(velocity[loc][0] * dNdx[loc][0] + velocity[loc][1] * dNdx[loc][1] + velocity[loc][2] * dNdx[loc][2]);
//		Nx[loc] += k_bar*(velocity[0] * dNdx[loc][0] + velocity[1] * dNdx[loc][1] + velocity[2] * dNdx[loc][2]);
//	}
//	detJ = dxdt.determinant();
//	detJ = 0.125*detJ;
//}
//
//void Laplace::ElementVelocity(vector<double>& Nx, const vector<array<double, 3>>& v_node, double v_tmp[3])
//{
//	int nen = v_node.size();
//	for (int j = 0; j < 3; j++)
//	{
//		v_tmp[j] = 0;
//	}
//	for (int i = 0; i < nen; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			v_tmp[j] += v_node[i][j] * Nx[i];
//		}
//	}
//}
//
//void Laplace::ElementMassMatrix(vector<double> &Weight, vector<double>& Nx, double detJ, vector<vector<double>>& EM)
//{
//	int i, j;
//	int nen = Nx.size();
//	for (i = 0; i<nen; i++)
//	{
//		for (j = 0; j<nen; j++)
//		{
//			EM[i][j] = Weight[i] * Nx[j] * detJ;
//		}
//	}
//}
//
//void Laplace::ElementConvectionMatrix(vector<double>& Nx, vector<array<double,3>>& dNdx, double v[3], double detJ, vector<vector<double>>& EC)
//{
//	int i, j;
//	int nen = Nx.size();
//	for (i = 0; i<nen; i++)
//	{
//		for (j = 0; j<nen; j++)
//		{
//			EC[i][j] = Nx[i] * (v[0] * dNdx[j][0] + v[1] * dNdx[j][1] + v[2] * dNdx[j][2])*detJ;
//		}
//	}
//}
//
//void Laplace::ElementStiffMatrix(vector<array<double,3>>& dNdx, double detJ, vector<vector<double>>& EK)
//{
//	int i, j;
//	int nen = dNdx.size();
//	for (i = 0; i<nen; i++)
//	{
//		for (j = 0; j<nen; j++)
//		{
//			EK[i][j] = (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
//		}
//	}
//
//}
//
//void Laplace::ElementLoadVector(vector<double>& Fx, double detJ, vector<double>& EF)
//{
//	int i;
//	int nen = Fx.size();
//	for (i = 0; i<nen; i++)
//	{
//		EF[i] = Fx[i] * detJ;
//	}
//}
//
//void Laplace::Tangent(const int nen, vector<vector<double>>& EM, vector<vector<double>>& EM_plus, vector<vector<double>>& EM_minus, vector<vector<double>>& EK, vector<vector<double>>& EC_plus, vector<vector<double>>& EC_minus, vector<vector<double>>& EMatrixSolve)
//{
//	int i, j;
//	for (i = 0; i < nen; i++)
//	{
//		for (j = 0; j < nen; j++)
//		{			
//			EMatrixSolve[i + 0][j + 0] += EM[i][j] + dt*DA*EK[i][j] + dt*(par[3] + par[4]) * EM[i][j];
//			EMatrixSolve[i + 0][j + nen] += -par[5] * dt * EM[i][j];
//			EMatrixSolve[i + nen][j + 0] += -par[3] * dt * EM_plus[i][j];
//			EMatrixSolve[i + nen][j + nen] += EM_plus[i][j] + dt*EC_plus[i][j] + dt*par[5] * EM_plus[i][j];
//
//			//EMatrixSolve[i + 0][j + 16] += -par[6] * dt * EM[i][j];
//			//EMatrixSolve[i + 16][j + 0] += -par[4] * dt * EM_minus[i][j];
//			//EMatrixSolve[i + 16][j + 16] += EM_minus[i][j] + dt*EC_minus[i][j] + dt*par[6] * EM_minus[i][j];
//		}
//	}
//}
//
//void Laplace::Residual(const int nen, vector<vector<double>>& EM, vector<vector<double>>& EM_plus, vector<vector<double>>& EM_minus, const vector<int>& IEN, vector<double> EVectorSolve)
//{
//	int i, j, k;
//	for (i = 0; i < nen; i++)
//	{
//		for (j = 0; j < nen; j++)
//		{
//			EVectorSolve[i] += EM[i][j] * CA[IEN[j]];
//			EVectorSolve[i + nen] += EM_plus[i][j] * N_plus[IEN[j]];
//			//EVectorSolve[i + nen*2] += EM_minus[i][j] * N_minus[IEN[j]];
//		}
//	}
//}

void Laplace::ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve)
{
	int j, k;
	int nen = EVectorSolve.size() / 2;
	for (j = 0; j < nen * 2; j++)
	{
		EVectorSolve[j] -= bc_value*EMatrixSolve[j][pt_num + variable_num * nen];
	}
	for (j = 0; j < nen * 2; j++)
	{
		EMatrixSolve[j][pt_num + variable_num *nen] = 0.0; EMatrixSolve[pt_num + variable_num * nen][j] = 0.0;
	}
	EMatrixSolve[pt_num + variable_num * nen][pt_num + variable_num * nen] = 1.0;
	EVectorSolve[pt_num + variable_num * nen] = bc_value;

}

void Laplace::Assembly(vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve, const vector<int>& IEN, vector<Triplet<double>>& trilist_MatrixSolve, vector<double>& VectorSolve)
{
	unsigned int i, j, A, B, p, q;
	int nen = IEN.size();
	for (i = 0; i<nen; i++)
	{
		A = IEN[i];
		for (j = 0; j<nen; j++)
		{
			B = IEN[j];
			if(EMatrixSolve[i + 0][j + 0]!=0.0)
			trilist_MatrixSolve.push_back(Triplet<double>(A, B, EMatrixSolve[i + 0][j + 0]));
			if (EMatrixSolve[i + 0][j + nen] != 0.0)
			trilist_MatrixSolve.push_back(Triplet<double>(A, B + CA.size(), EMatrixSolve[i + 0][j + nen]));
			if (EMatrixSolve[i + nen][j + 0] != 0.0)
			trilist_MatrixSolve.push_back(Triplet<double>(A + CA.size(), B, EMatrixSolve[i + nen][j + 0]));
			if (EMatrixSolve[i + nen][j + nen] != 0.0)
			trilist_MatrixSolve.push_back(Triplet<double>(A + CA.size(), B + CA.size(), EMatrixSolve[i + nen][j + nen]));
			//trilist_MatrixSolve.push_back(Triplet<double>(A, B + CA.size() * 2, EMatrixSolve[i + 0][j + 16]));
			//trilist_MatrixSolve.push_back(Triplet<double>(A + CA.size() * 2, B, EMatrixSolve[i + 16][j + 0]));
			//trilist_MatrixSolve.push_back(Triplet<double>(A + CA.size() * 2, B + CA.size() * 2, EMatrixSolve[i + 16][j + 16]));

		}
	}
	for (i = 0; i < nen; i++)
	{
		A = IEN[i];
		VectorSolve[A] += EVectorSolve[i];
		VectorSolve[A + CA.size()] += EVectorSolve[i + nen];
		//VectorSolve[A + CA.size() * 2] += EVectorSolve[i + 16];

	}
}

void Laplace::BuildLinearSystemIGA(const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<array<double, 3>> &cpts, const vector<int>& pid_loc, const vector<int>& label, const vector<array<double, 3>> velocity_node, const double Vplus, const double Vminus, SparseMatrix<double, RowMajor>& MatrixSolve, vector<double>& VectorSolve)
{

#pragma omp declare reduction (merge : vector<Triplet<double>>: omp_out.insert(omp_out.end(), omp_in.begin(),omp_in.end()))
	vector<Triplet<double>> trilist_MatrixSolve;
	trilist_MatrixSolve.clear();

#pragma omp parallel for reduction(merge: trilist_MatrixSolve)
	for (int e = 0; e<bzmesh.size(); e++)
	{
		//cout << "Element: " << e << endl;
		double detJ;
		double tau[2];
		double dudx[3][3];
		int nen = bzmesh[e].IEN.size();

		vector<vector<double>> EM, EK, EC_plus, EC_minus, EM_plus, EM_minus;
		vector<vector<double>> EMatrixSolve;
		vector<double>EVectorSolve;
		vector<double>Nx, Npx, Nmx;
		vector<array<double, 3>> dNdx;

		EM.clear();		EK.clear();		EC_plus.clear();	EC_minus.clear();	EM_plus.clear();	EM_minus.clear();
		EMatrixSolve.clear();	EVectorSolve.clear();
		Nx.clear(); Npx.clear(); Nmx.clear();
		dNdx.clear();
		
		EM.resize(nen);		EK.resize(nen);		EC_plus.resize(nen);	EC_minus.resize(nen);	EM_plus.resize(nen);	EM_minus.resize(nen);
		EMatrixSolve.resize(nen * 2);	EVectorSolve.resize(nen*2);
		Nx.resize(nen); Npx.resize(nen); Nmx.resize(nen); 
		dNdx.resize(nen);
		
		for (int i = 0; i < nen; i++)
		{
			EM[i].resize(nen);		EK[i].resize(nen);		EC_plus[i].resize(nen);		EC_minus[i].resize(nen);	EM_plus[i].resize(nen);	EM_minus[i].resize(nen);
			EMatrixSolve[i].resize(nen * 2, 0.0);
			EMatrixSolve[i+nen].resize(nen * 2, 0.0);
		}

		for (int i = 0; i < nen * 2; i++)
		{
			for (int j = 0; j < nen * 2; j++)
			{
				EMatrixSolve[i][j] = 0.0;
			}
			EVectorSolve[i] = 0.0;
		}

		double vplus_tmp[3], vminus_tmp[3];
		double CA_tmp, Nplus_tmp, Nminus_tmp;

		vector<double> CA_cpt, Nplus_cpt, Nminus_cpt;
		CA_cpt.resize(nen); Nplus_cpt.resize(nen); Nminus_cpt.resize(nen);
		vector<array<double, 3>> vplus_cpt, vminus_cpt;
		vplus_cpt.resize(nen); vminus_cpt.resize(nen);

		for (int i = 0; i < nen; i++)
		{
			CA_cpt[i] = CA[bzmesh[e].IEN[i]];
			Nplus_cpt[i] = N_plus[bzmesh[e].IEN[i]];
			Nminus_cpt[i] = N_minus[bzmesh[e].IEN[i]];
			for (int j = 0; j < dim; j++)
			{
				vplus_cpt[i][j] = Vplus*velocity_node[bzmesh[e].IEN[i]][j];
				vminus_cpt[i][j] = Vminus*velocity_node[bzmesh[e].IEN[i]][j];
			}
		}

		for (int i = 0; i<Gpt.size(); i++)
		{
			for (int j = 0; j<Gpt.size(); j++)
			{
				for (int k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, bzmesh[e].cmat, Nx, dNdx, dudx, detJ);
					ElementVelocity(Nx, vplus_cpt, vplus_tmp);
					ElementVelocity(Nx, vminus_cpt, vminus_tmp);
					ElementValue(Nx, CA_cpt, CA_tmp);
					ElementValue(Nx, Nplus_cpt, Nplus_tmp);
					ElementValue(Nx, Nminus_cpt, Nminus_tmp);
					SUPGcoefficientIGA(par[5], vplus_tmp, dNdx, tau[0]);
					SUPGcoefficientIGA(par[6], vminus_tmp, dNdx, tau[1]);
					WeightingFunction(vplus_tmp, tau[0], Nx, dNdx, Npx);
					WeightingFunction(vminus_tmp, tau[1], Nx, dNdx, Nmx);
					detJ = wght[i] * wght[j] * wght[k] * detJ;

					ElementMassMatrix(Nx, Nx, detJ, EM);
					ElementMassMatrix(Npx, Nx, detJ, EM_plus);
					ElementMassMatrix(Nmx, Nx, detJ, EM_minus);
					ElementStiffMatrix(dNdx, detJ, EK);
					ElementConvectionMatrix(Npx, dNdx, vplus_tmp, detJ, EC_plus);
					ElementConvectionMatrix(Nmx, dNdx, vminus_tmp, detJ, EC_minus);
					Tangent(nen, EM, EM_plus, EM_minus, EK, EC_plus, EC_minus, EMatrixSolve);
					Residual(nen, CA_tmp, Nplus_tmp, Nminus_tmp, Nx, Npx, Nmx, detJ, EVectorSolve);
				}
			}
		}

		//Apply BCs
		double CA_bc, Nplus_bc, Nminus_bc;
		for (int i = 0; i < nen; i++)
		{
			int A = bzmesh[e].IEN[i];
			////Cylinder BCs			
			//if (sqrt(pow(cpts[A][1], 2) + pow(cpts[A][2], 2)) > 0.48)
			//{
			//	CA_bc = 0;
			//	ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//	//Nminus_bc = 0;
			//	//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//}
			//else if (label[A] == 1 || label[A] == 2)

			//Cylinder photoactivation BCs
			//if (label[A] == 1)
			//{
			//	CA_bc = 0.0;
			//	ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	Nplus_bc = 0.0;
			//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//}

			////Bifurcation BCs
			////if (pid_loc[A] != -1)
			//{
			//	//double R0 = sqrt(pow(cpts[A][0] - 0.0, 2) + pow(cpts[A][1] - 0.0, 2) + pow(cpts[A][2] - 0.0, 2));
			//	//double R1 = sqrt(pow(cpts[A][0] - 7.33177, 2) + pow(cpts[A][1] - (-6.29269), 2) + pow(cpts[A][2] - (-5.655), 2));
			//	//double R2 = sqrt(pow(cpts[A][0] - 4.85463, 2) + pow(cpts[A][1] - (-5.84302), 2) + pow(cpts[A][2] - 1.701, 2));
			//
			//	if (label[A] == 1)
			//	{
			//		//if (R0 > 0.7)//inlet wall
			//		//{
			//		//	CA_bc = 0.0;
			//		//	ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		//	Nplus_bc = 0.0;
			//		//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//	//Nminus_bc = 0;
			//		//	//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		//}
			//		//else//inlet
			//		{
			//			CA_bc = 1.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//			Nplus_bc = 2.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//
			//	}
			//	if (label[A] == 0) //wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	//if (label[A] == 2 && R1>0.45) //outlet1 wall
			//	//{
			//	//	CA_bc = 0;
			//	//	ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	//	Nplus_bc = 0;
			//	//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//	//	//Nminus_bc = 0;
			//	//	//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	//}
			//	//if (label[A] == 3 && R2>0.47) //outlet2 wall
			//	//{
			//	//	CA_bc = 0;
			//	//	ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//	//	Nplus_bc = 0;
			//	//	ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//	//	//Nminus_bc = 0;
			//	//	//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	//}
			//
			//}

			////Bifurcation photoactivation BCs
			if ((A<201))
			{
				CA_bc = 0.0;
				ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
				Nplus_bc = 0.0;
				ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			}
			
			////3 Bifurcation BCs
			//if (pid_loc[A] != -1)
			//{
			//	double R0 = sqrt(pow(cpts[A][0] - 0.0, 2)     + pow(cpts[A][1] - 0.0, 2)     + pow(cpts[A][2] - 0.0, 2));
			//	double R1 = sqrt(pow(cpts[A][0] - 63.6235, 2) + pow(cpts[A][1] - 14.7335, 2) + pow(cpts[A][2] - 8.25363, 2));
			//	double R2 = sqrt(pow(cpts[A][0] - 46.7595, 2) + pow(cpts[A][1] - 22.4332, 2) + pow(cpts[A][2] - 7.46334, 2));
			//	double R3 = sqrt(pow(cpts[A][0] - 35.7126, 2) + pow(cpts[A][1] - 36.7129, 2) + pow(cpts[A][2] - 2.66824, 2));
			//	double R4 = sqrt(pow(cpts[A][0] - 28.7469, 2) + pow(cpts[A][1] - 43.6465, 2) + pow(cpts[A][2] - 2.81018, 2));
			//	
			//	if (label[A] == 1)
			//	{
			//		if (R0 > 2.0)//inlet wall
			//		{
			//			CA_bc = 0.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//			Nplus_bc = 0.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//		else//inlet
			//		{
			//			CA_bc = 1.0;
			//			ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//			Nplus_bc = 2.0;
			//			ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//			//Nminus_bc = 0;
			//			//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//		}
			//
			//	}
			//	if (label[A] < 1) //wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 2 && R1>0.78) //outlet1 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 3 && R2>0.61) //outlet2 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 4 && R3>0.64) //outlet3 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//	if (label[A] == 5 && R4>0.71) //outlet4 wall
			//	{
			//		CA_bc = 0;
			//		ApplyBoundaryCondition(CA_bc, i, 0, EMatrixSolve, EVectorSolve);
			//		Nplus_bc = 0;
			//		ApplyBoundaryCondition(Nplus_bc, i, 1, EMatrixSolve, EVectorSolve);
			//		//Nminus_bc = 0;
			//		//ApplyBoundaryCondition(Nminus_bc, i, 2, EMatrixSolve, EVectorSolve);
			//	}
			//
			//}
		}
		Assembly(EMatrixSolve, EVectorSolve, bzmesh[e].IEN, trilist_MatrixSolve, VectorSolve);
	}
	//cout << "Assemble Global Matrix" << endl;
	if (judge == 0)
		MatrixSolve.setFromTriplets(trilist_MatrixSolve.begin(), trilist_MatrixSolve.end());
}

//void Laplace::Solver(const vector<array<double, 3>>& pts, const vector<int>& label, const vector<int>& pid_loc,double time)
//{
//	cout << "Diffusion...\n";
//	VectorXd tmp_Nplus(N_plus.size()), tmp_Nminus(N_minus.size()), tmp_CA(CA.size());
//	VectorXd GF_mix(CA.size() * 2);
//	for (uint i = 0; i < CA.size(); i++)
//	{
//		tmp_CA[i] = CA[i];
//	}
//	for (uint i = 0; i < N_plus.size(); i++)
//	{
//		tmp_Nplus[i] = N_plus[i];
//	}
//	for (uint i = 0; i < N_minus.size(); i++)
//	{
//		tmp_Nminus[i] = N_minus[i];
//	}
//	VectorXd GF_CA = GM*tmp_CA;
//	VectorXd GF_plus = GM_plus*tmp_Nplus;
//	VectorXd GF_minus = GM_minus*tmp_Nminus;
//
//	
//	SparseMatrix<double> GM1(CA.size(), CA.size());
//	SparseMatrix<double> GM_plus1(CA.size(), CA.size());
//	SparseMatrix<double> GM_minus1(CA.size(), CA.size());
//	SparseMatrix<double> GMtest(CA.size(), CA.size());
//	
//	GM1 = GM + dt* DA* GK + dt* (par[3] + par[4]) * GM;	
//	GM_plus1 = GM_plus + dt* GC_plus + dt* par[5] * GM_plus;
//	GM_minus1 = GM_minus + dt* GC_minus + dt* par[6] * GM_minus;
//	GMtest = GM + dt* DA* GK + dt* GC_plus + dt* (par[3] + par[4]) * GM;//convection-diffusion test
//
//
//
//	//Boundary Conditions
//	for (uint i = 0; i < label.size(); i++)
//	{
//		if (label[i] == 1)
//		{
//			GF_CA[i] = 1.*1.e20*GM1.coeff(i, i);
//			GF_plus[i] = 0.2 *1.e20* GM_plus1.coeff(i, i);
//			//GF_CA[i] = 1.*1.e11*GMtest.coeff(i, i);
//		}
//		else if (label[i] > 1)
//		{
//			//GF_CA[i] = 0.3*1.e11*GM1.coeff(i, i);
//			//GF_CA[i] = 0.3*1.e11*GMtest.coeff(i, i);
//			//GF_minus[i] = 0. * 1.e11* GM_minus1.coeffRef(i, i);
//		}
//		if (sqrt(pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]) > 0.48)
//		{
//			GF_CA[i] = 0;
//			GF_plus[i] = 0;
//			//GF_minus[i] = 0;
//		}
//	}
//	
//	//string fn55;
//	//fn55 = "../io/matrix_examine/GF_CA.txt";
//	//ofstream fout55;
//	//fout55.open(fn55);
//	//for(uint i=0;i<label.size();i++)
//	//	fout55 << GF_CA[i]<<"\n";
//
//	for (uint i = 0; i < CA.size(); i++)
//	{
//		GF_mix[i] = GF_CA[i];
//		GF_mix[i + CA.size()] = GF_plus[i];
//		//GF_mix[i + CA.size() * 2] = GF_minus[i];
//	}
//	cout << "Solving...\n";
//
//
//	VectorXd sol_mix = solver_mix.solve(GF_mix);
//	for (uint i = 0; i < CA.size(); i++)
//	{
//		CA[i] = sol_mix[i];
//		N_plus[i] = sol_mix[i + CA.size()];
//		//N_minus[i] = sol_mix[i + CA.size() * 2];
//	}
//	cout << "Done diffusion!\n";
//}

void Laplace::ConcentrationCal(double u, double v, double w, const Element3D& bzel, double& disp, double& detJ)
{
	double Nx[64];
	double dNdx[64][3];
	double dudx[3][3];
	BasisFunction(u, v, w, bzel.pts, Nx, dNdx, dudx, detJ);
	double uloc[64];
	unsigned int i, j;
	for (i = 0; i<64; i++)
	{
		uloc[i] = 0.;
		for (j = 0; j<bzel.IEN.size(); j++)
		{
			uloc[i] += bzel.cmat[j][i] * (CA[bzel.IEN[j]]+N_plus[bzel.IEN[j]]);
		}
	}
	//displacement
	disp = 0.;
	for (i = 0; i<64; i++)
	{
		disp += Nx[i] * uloc[i];
	}
}

void Laplace::ConcentrationCal_Coupling_Bezier(double u, double v, double w, const Element3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ)
{
	double dUdx[3][3];
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 3>> dNdx(bzel.IEN.size());
	bzel.Para2Phys(u, v, w, pt);
	BasisFunction(u, v, w, bzel.pts, bzel.cmat, Nx, dNdx, dUdx, detJ);
	disp = 0.;
	dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)
	{
		disp += Nx[i] * (CA[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[0] += dNdx[i][0] * (CA[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[1] += dNdx[i][1] * (CA[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
		dudx[2] += dNdx[i][2] * (CA[bzel.IEN[i]] + N_plus[bzel.IEN[i]]);
	}
}

void Laplace::VisualizeVTK(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn)
{
	string fname;
	ofstream fout;
	unsigned int i;

	//fname = fn + "_N0.vtk";
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
	//	fout<<"POINT_DATA "<<CA.size()<<"\nSCALARS N0 float 1\nLOOKUP_TABLE default\n";
	//	for(uint i=0; i<CA.size(); i++)
	//	{
	//		fout<<CA[i]<<"\n";
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}
	//
	//fname = fn + "_Nplus.vtk";
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

	//fname = fn + "_Nminus.vtk";
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
	//	fout << "POINT_DATA " << N_minus.size() << "\nSCALARS N- float 1\nLOOKUP_TABLE default\n";
	//	for (uint i = 0; i<N_minus.size(); i++)
	//	{
	//		fout << N_minus[i] << "\n";
	//	}
	//	fout.close();
	//}
	//else
	//{
	//	cout << "Cannot open " << fname << "!\n";
	//}

	fname = fn + "_allparticle.vtk";
	fout.open(fname.c_str());
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
		fout << "POINT_DATA " << N_plus.size() << "\nSCALARS AllParticles float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<N_plus.size(); i++)
		{
			fout << CA[i]+N_plus[i]+N_minus[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Laplace::VisualizeVTK_1(const vector<Element3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	//vector<array<double, 3>> sdisp_err;
	//vector<array<double, 3>> sse;
	//vector<array<double, 3>> sss;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
							  //int ecount(0);
							  //vector<double> su;
							  //for (int i = 0; i<ns; i++)
							  //{
							  //	su[i] = double(i) / (double(ns) - 1.);
							  //}

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		//double errtmp;
		//ElementError(bzmesh[e],errtmp);
		//errL2.push_back(sqrt(errtmp));

		//if (bzmesh[e].type == 1)
		{
			int ns(2);
			if (bzmesh[e].type == 1) ns = 5;
			vector<double> su(ns);
			for (int i = 0; i<ns; i++)
			{
				su[i] = double(i) / (double(ns) - 1.);
			}

			int loc(0);
			int pstart = spt.size();
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						double pt1[3];
						double disp, detmp;
						bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
						ConcentrationCal(su[c], su[b], su[a], bzmesh[e], disp, detmp);
						array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
						spt.push_back(pt);
						sdisp.push_back(disp);
					}
				}
			}
			int nns[2] = { ns*ns*ns, ns*ns };
			for (int a = 0; a<ns - 1; a++)
			{
				for (int b = 0; b<ns - 1; b++)
				{
					for (int c = 0; c < ns - 1; c++)
					{
						array<int, 8> el;
						el[0] = pstart + a*nns[1] + b*ns + c;
						el[1] = pstart + a*nns[1] + b*ns + c + 1;
						el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
						el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
						el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
						el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
						el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
						el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
						sele.push_back(el);
					}
				}
			}
		}
	}

	string fname = fn + "_allparticle.vtk";
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
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Laplace::VisualizeVTK_2(const vector<Element3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	double detJ;

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int ns(3); //sample points number
		//if (bzmesh[e].type == 1) ns = 5;
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		int pstart = spt.size();
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3], dudx[3];
					double disp;
					ConcentrationCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh[e], pt1, disp, dudx, detJ);
					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
					spt.push_back(pt);
					sdisp.push_back(disp);
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = pstart + a*nns[1] + b*ns + c;
					el[1] = pstart + a*nns[1] + b*ns + c + 1;
					el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
					el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
					el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
					el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
	}

	string fname = fn + "_allparticle.vtk";
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
		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS AllParticle float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Laplace::Solver(const vector<array<double, 3>>& cpts, const vector<array<double, 3>> velocity_node, const vector<Element3D>& bzmesh, const vector<Element3D>& tmesh, const vector<int>& label, const vector<int>& pid_loc)
{
	cout << "Diffusion...\n";
	GaussInfo(4);
	VectorXd tmp_solution(2 * cpts.size()), tmp_VectorSolve(2 * cpts.size());


	cout << "Building linear system...\n";
	VectorSolve.clear();
	VectorSolve.resize(2 * cpts.size(), 0.0);
	MatrixSolve.setZero();
	BuildLinearSystem(bzmesh, tmesh, cpts, pid_loc, label, velocity_node, par[1], par[2], MatrixSolve, VectorSolve);


	for (uint i = 0; i < VectorSolve.size(); i++)
	{
		tmp_VectorSolve[i] = VectorSolve[i];
	}


	//cout << "Solving...\n";

	tmp_solution = solver_MIX.solve(tmp_VectorSolve);

	for (uint i = 0; i < CA.size(); i++)
	{
		CA[i] = tmp_solution[i];
		N_plus[i] = tmp_solution[i + CA.size()];
		//N_minus[i] = tmp_solution[i + CA.size() * 2];
	}


	//cout << "Done diffusion!\n";

}

//void Laplace::PardisoSolver(SparseMatrix<double, RowMajor>& GK, vector<double>& GF, vector<double>& solution)
//{
//	/* Matrix data. */
//	MKL_INT n = GK.rows();
//	vector<MKL_INT> JA;
//	vector<double> A;
//	MKL_INT n_nonzero = GK.nonZeros();
//	MKL_INT *ia = new MKL_INT[n + 1];
//	//MKL_INT *ja = new MKL_INT[n_nonzero];
//	//double *a = new double[n_nonzero];
//	int ii, jj, index(0), flag_first(0);
//	for (ii = 0; ii < GK.outerSize(); ii++)
//	{
//		for (SparseMatrix<double, RowMajor>::InnerIterator it(GK, ii); it; ++it)
//		{
//			if (flag_first == 0)
//			{
//				ia[ii] = index + 1;
//				flag_first = 1;
//			}
//			JA.push_back(it.col() + 1);
//			A.push_back(it.value());
//			index = index + 1;
//		}
//		flag_first = 0;
//	}
//	n_nonzero = A.size();
//	ia[n] = n_nonzero + 1;
//	MKL_INT *ja = new MKL_INT[n_nonzero];
//	double *a = new double[n_nonzero];
//
//	for (int ii = 0; ii < n_nonzero; ii++)
//	{
//		ja[ii] = JA[ii];
//		a[ii] = A[ii];
//	}
//	
//	MKL_INT mtype = 11;       /* Real unsymmetric matrix */
//							  /* RHS and solution vectors. */
//	double *b = new double[n];
//	double *x = new double[n];
//	double *bs = new double[n];
//	double res, res0;
//	MKL_INT nrhs = 1;     /* Number of right hand sides. */
//						  /* Internal solver memory pointer pt, */
//						  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
//						  /* or void *pt[64] should be OK on both architectures */
//	void *pt[64];
//	/* Pardiso control parameters. */
//	MKL_INT iparm[64];
//	MKL_INT maxfct, mnum, phase, error, msglvl;
//	/* Auxiliary variables. */
//	MKL_INT i, j;
//	double ddum;          /* Double dummy */
//	MKL_INT idum;         /* Integer dummy. */
//	char *uplo;
//	/* -------------------------------------------------------------------- */
//	/* .. Setup Pardiso control parameters. */
//	/* -------------------------------------------------------------------- */
//	for (i = 0; i < 64; i++)
//	{
//		iparm[i] = 0;
//	}
//	iparm[0] = 1;         /* No solver default */
//	iparm[1] = 2;         /* Fill-in reordering from METIS */
//	iparm[3] = 0;         /* No iterative-direct algorithm */
//	iparm[4] = 0;         /* No user fill-in reducing permutation */
//	iparm[5] = 0;         /* Write solution into x */
//	iparm[6] = 0;         /* Not in use */
//	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
//	iparm[8] = 0;         /* Not in use */
//	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
//	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
//	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
//	iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
//	iparm[13] = 0;        /* Output: Number of perturbed pivots */
//	iparm[14] = 0;        /* Not in use */
//	iparm[15] = 0;        /* Not in use */
//	iparm[16] = 0;        /* Not in use */
//	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
//	iparm[18] = -1;       /* Output: Mflops for LU factorization */
//	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
//	iparm[26] = 1;
//	maxfct = 1;           /* Maximum number of numerical factorizations. */
//	mnum = 1;         /* Which factorization to use. */
//	msglvl = 1;           /* Print statistical information in file */
//	error = 0;            /* Initialize error flag */
//						  /* -------------------------------------------------------------------- */
//						  /* .. Initialize the internal solver memory pointer. This is only */
//						  /* necessary for the FIRST call of the PARDISO solver. */
//						  /* -------------------------------------------------------------------- */
//	for (i = 0; i < 64; i++)
//	{
//		pt[i] = 0;
//	}
//	
//	/* -------------------------------------------------------------------- */
//	/* .. Reordering and Symbolic Factorization. This step also allocates */
//	/* all memory that is necessary for the factorization. */
//	/* -------------------------------------------------------------------- */
//	phase = 11;
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0)
//	{
//		printf("\nERROR during symbolic factorization: %d", error);
//		getchar();
//		exit(1);
//	}
//	printf("\nReordering completed ... ");
//	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//	/* -------------------------------------------------------------------- */
//	/* .. Numerical factorization. */
//	/* -------------------------------------------------------------------- */
//	phase = 22;
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0)
//	{
//		printf("\nERROR during numerical factorization: %d", error);
//		getchar();
//		exit(2);
//	}
//	printf("\nFactorization completed ... "); 
//	
//	/* -------------------------------------------------------------------- */
//	/* .. Back substitution and iterative refinement. */
//	/* -------------------------------------------------------------------- */
//	phase = 33;
//	/* Set right hand side to one. */
//	for (i = 0; i < n; i++)
//	{
//		b[i] = GF[i];
//	}
//	//uplo = "non-transposed";
//
//
//	printf("\n\nSolving %s system...\n", uplo);
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
//	if (error != 0)
//	{
//		printf("\nERROR during solution: %d", error);
//		getchar();
//		exit(3);
//	}
//
//	for (j = 0; j < n; j++)
//	{
//		solution[j] = x[j];
//	}
//	//printf("\nThe solution of the system is: ");
//	//for (j = 0; j < n; j++)
//	//{
//	//	printf("\n x [%d] = % f", j, x[j]);
//	//}
//	printf("\n");
//	// Compute residual
//	mkl_dcsrgemv(uplo, &n, a, ia, ja, x, bs);
//	res = 0.0;
//	res0 = 0.0;
//	for (j = 1; j <= n; j++)
//	{
//		res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
//		res0 += b[j - 1] * b[j - 1];
//	}
//	res = sqrt(res) / sqrt(res0);
//	printf("\nRelative residual = %e", res);
//	// Check residual
//	if (res > 1e-10)
//	{
//		printf("Error: residual is too high!\n");
//		exit(10 + i);
//	}
//
//
//	/* -------------------------------------------------------------------- */
//	/* .. Termination and release of memory. */
//	/* -------------------------------------------------------------------- */
//	phase = -1;           /* Release internal memory. */
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, &ddum, ia, ja, &idum, &nrhs,
//		iparm, &msglvl, &ddum, &ddum, &error);
//	printf("\n");
//}

void Laplace::Run(const vector<array<double, 3>>& cpts, const vector<array<double, 3>> velocity_node, const vector<Element3D> &tmesh, vector<Element3D>& bzmesh, const vector<int>& label, const vector<int>& pid_loc)
{
	cout << "Diffusion...\n";
	GaussInfo(4);
	VectorXd tmp_solution(2 * cpts.size()), tmp_VectorSolve(2 * cpts.size());
	//VectorXd tmp_VectorSolve(2 * cpts.size());
	//vector<double> tmp_solution(2 * cpts.size());

	cout<<"Building linear system...\n";
	VectorSolve.clear();
	VectorSolve.resize(2 * cpts.size(), 0.0);
	MatrixSolve.setZero();
	//BuildLinearSystem(bzmesh, tmesh, cpts, pid_loc, label, velocity_node, par[1], par[2], MatrixSolve, VectorSolve);
	BuildLinearSystemIGA(bzmesh, tmesh, cpts, pid_loc, label, velocity_node, par[1], par[2], MatrixSolve, VectorSolve);
	MatrixSolve.makeCompressed();

	for (uint i = 0; i < VectorSolve.size(); i++)
	{
		tmp_VectorSolve[i] = VectorSolve[i];
	}
	
	cout << "Solving...\n";
	PardisoLU<SparseMatrix<double>> solver0;
	solver0.compute(MatrixSolve);
	tmp_solution = solver0.solve(tmp_VectorSolve);

	if (judge == 0)
	{
		solver_MIX.compute(MatrixSolve);
		judge = 1;
	}
	tmp_solution = solver_MIX.solve(tmp_VectorSolve);

	//PardisoSolver(MatrixSolve, VectorSolve, tmp_solution);

	for (uint i = 0; i < CA.size(); i++)
	{
		CA[i] = tmp_solution[i];
		N_plus[i] = tmp_solution[i + CA.size()];
		//N_minus[i] = tmp_solution[i + CA.size() * 2];
	}
	cout << "Done diffusion!\n";

	//for (int e = 0; e < bzmesh.size(); e++)
	//{
	//	Tmesh2Bzmesh_value(bzmesh[e], CA, bzmesh[e].CA_bz);
	//	Tmesh2Bzmesh_value(bzmesh[e], N_plus, bzmesh[e].Nplus_bz);
	//	Tmesh2Bzmesh_value(bzmesh[e], N_minus, bzmesh[e].Nminus_bz);
	//}

	//string fn44;
	//fn44 = "../io/matrix_examine/GM1.txt";
	//ofstream fout44;
	//fout44.open(fn44);
	//fout44 << MatrixSolve1;
}


