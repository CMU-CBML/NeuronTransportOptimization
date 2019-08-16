#ifndef DIFF_REACT_H
#define DIFF_REACT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
#include <cmath>
using namespace std;



//problem setting

void SetVariables(vector<double>& var);

void ReadMesh(string fn, vector<array<double, 3>>& pts, vector<Element3D>& mesh, vector<array<double, 3>>& pts_b, vector<Element2D>& mesh_b, vector<int>& pid_loc);

void ReadMeshLabel(string fn, vector<array<double, 3>>& pts, vector<int>& label, vector<Element3D>& mesh, vector<array<double, 3>>& pts_b, vector<Element2D>& mesh_b, vector<int>& pid_loc);

void ReadBezierElement(string fn, vector<Element3D>& mesh);

void ReadVelocityField(string fn, vector<Element3D>& mesh);

void ReadVelocityFieldNode(string fn, vector<array<double, 3>>& pts, vector<array<double, 3>>& velocity);

void SetInitialCondition(int ndof, int ndof_b, vector<double>& CA0, vector<double>& NX0, vector<double>& NB0, vector<array<double, 3>>& pts, const vector<int>& label, const vector<int>& pid_loc);

void SetTestBC(double& CAi, double& CAs, double& length);

void SetTempData(int ndof, int ndof_b, vector<double>& CA, vector<double>& NX, vector<double>& NB);

double MaxDifference(vector<double> a, vector<double> b);

void DataTrans2React(int ndof_b, const vector<int>& pid_loc, const vector<double>& CA, vector<double>& CA_b);

void DataTrans2DiffuseBv(int ndof, const vector<int>& pid_loc, const vector<double>& Bv, vector<double>& Bv_all);

void DataTrans2Diffuse(int ndof, const vector<int>& pid_loc, const vector<double>& CA_b, vector<double>& CA);

#endif