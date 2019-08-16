#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include "mkl.h"
#include <Eigen/Dense>
#include <Eigen/PardisoSupport>
using namespace std;
using namespace Eigen;


const int bzpt_num = 64;
const int degree = 3;
const int dim = 3;
/////////////////////////////////

class Element2D
{
public:
	vector<int> cnct;
	vector<int> IEN;
	vector<vector<double>> cmat;
	vector<array<double, 3>> pts;//tmp

	Element2D();
};

class Element3D
{
public:
	int degree;
	int order;
	int nbf;
	int type;//0 for interior and 1 for boundary, for visualization purpose
	int bzflag;//0 for spline element, 1 for Bezier element

	vector<int> IEN;
	vector<int> IENb;
	vector<array<double, 64>> cmat;
	vector<array<double, 3>> pts;//tmp
	vector<array<double, 3>> vplus_bz;
	vector<array<double, 3>> vminus_bz;
	double velocity[3];
	double CA_bz[bzpt_num];
	double Nplus_bz[bzpt_num];
	double Nminus_bz[bzpt_num];

	Element3D(int p = 3);
	void BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const;
	void Basis(double u, double v, double w, vector<double>& Nt, vector<array<double, 3>>& dNdt) const;
	void Para2Phys(double u, double v, double w, double pt[3]) const;
	
};

class ShapeFunction
{
public:
	const int Num_Gpt = 64;
	vector<array<double, 64>> Nx;
	vector<array<array<double, 3>, 64>> dNdx;

	ShapeFunction();
};

//mesh format conversion

void Raw2Vtk_hex(string fn);

void Rawn2Vtk_hex(string fn);

#endif