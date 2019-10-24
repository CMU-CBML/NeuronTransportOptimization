#include "BasicDataStructure.h"

Vertex2D::Vertex2D()
{
	coor[0] = 0.;	coor[1] = 0.;
	label = 0;
}

Element2D::Element2D(int p)
{
	degree = p;
	order = p + 1;
	nbf = order*order*order;
	IEN.resize(4);
	pts.resize(4);
	for (int i = 0; i < 4; i++)
	{
		IEN[i] = 0;
		pts[i][0] = 0.; pts[i][1] = 0.;
	}
}

void Element2D::BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const
{
	if (degree == 3)
	{
		double Nu0[4] = { (1. - u)*(1. - u)*(1. - u),3.*(1. - u)*(1. - u)*u,3.*(1. - u)*u*u,u*u*u };
		double dNdu0[4] = { -3.*(1. - u)*(1. - u),3. - 12.*u + 9.*u*u,3.*(2. - 3.*u)*u,3.*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
	else if (degree == 4)
	{
		double Nu0[5] = { (1. - u)*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - u)*u,6.*(1. - u)*(1. - u)*u*u,4.*(1. - u)*u*u*u,u*u*u*u };
		double dNdu0[5] = { -4.*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - 4.*u),12.*u*(1. - 3.*u + 2.*u*u),4.*(3. - 4.*u)*u*u,4.*u*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
}

void Element2D::Basis(double u, double v, vector<double>& Nt, vector<array<double, 2>>& dNdt) const
{
	vector<double> Nu, Nv, dNdu, dNdv;
	BezierPolyn(u, Nu, dNdu);
	BezierPolyn(v, Nv, dNdv);
	Nt.resize(nbf);
	dNdt.resize(nbf);
	int i, j, loc(0);
	for (j = 0; j<order; j++)
	{
		for (i = 0; i<order; i++)
		{
			Nt[loc] = Nu[i] * Nv[j];
			dNdt[loc][0] = dNdu[i] * Nv[j];
			dNdt[loc][1] = Nu[i] * dNdv[j];
			loc++;
		}
	}
}

void Element2D::Para2Phys(double u, double v, double pt[2]) const
{
	vector<double> Nt;
	vector<array<double, 2>> dNdt;
	Basis(u, v, Nt, dNdt);
	pt[0] = 0.; pt[1] = 0.;
	for (int i = 0; i<nbf; i++)
	{
		pt[0] += pts[i][0] * Nt[i];
		pt[1] += pts[i][1] * Nt[i];
	}
}


Vertex3D::Vertex3D()
{
	coor[0] = 0.;	coor[1] = 0.;	coor[2] = 0.;
	label = 0;
}

Element3D::Element3D(int p)
{
	degree = p;
	order = p + 1;
	nbf = order*order*order;
	IEN.resize(8);
	pts.resize(8);
	for (int i = 0; i < 8; i++)
	{
		IEN[i] = 0;
		pts[i][0] = 0.; pts[i][1] = 0.; pts[i][2] = 0.;
	}
}

void Element3D::BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const
{
	if (degree == 3)
	{
		double Nu0[4] = { (1. - u)*(1. - u)*(1. - u),3.*(1. - u)*(1. - u)*u,3.*(1. - u)*u*u,u*u*u };
		double dNdu0[4] = { -3.*(1. - u)*(1. - u),3. - 12.*u + 9.*u*u,3.*(2. - 3.*u)*u,3.*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
	else if (degree == 4)
	{
		double Nu0[5] = { (1. - u)*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - u)*u,6.*(1. - u)*(1. - u)*u*u,4.*(1. - u)*u*u*u,u*u*u*u };
		double dNdu0[5] = { -4.*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - 4.*u),12.*u*(1. - 3.*u + 2.*u*u),4.*(3. - 4.*u)*u*u,4.*u*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
}

void Element3D::Basis(double u, double v, double w, vector<double>& Nt, vector<array<double, 3>>& dNdt) const
{
	vector<double> Nu, Nv, Nw, dNdu, dNdv, dNdw;
	BezierPolyn(u, Nu, dNdu);
	BezierPolyn(v, Nv, dNdv);
	BezierPolyn(w, Nw, dNdw);
	Nt.resize(nbf);
	dNdt.resize(nbf);
	int i, j, k, loc(0);
	for (k = 0; k<order; k++)
	{
		for (j = 0; j<order; j++)
		{
			for (i = 0; i<order; i++)
			{
				Nt[loc] = Nu[i] * Nv[j] * Nw[k];
				dNdt[loc][0] = dNdu[i] * Nv[j] * Nw[k];
				dNdt[loc][1] = Nu[i] * dNdv[j] * Nw[k];
				dNdt[loc][2] = Nu[i] * Nv[j] * dNdw[k];
				loc++;
			}
		}
	}
}

void Element3D::Para2Phys(double u, double v, double w, double pt[3]) const
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	Basis(u, v, w, Nt, dNdt);
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (int i = 0; i<nbf; i++)
	{
		pt[0] += pts[i][0] * Nt[i];
		pt[1] += pts[i][1] * Nt[i];
		pt[2] += pts[i][2] * Nt[i];
	}
}
