#ifndef UTILS_H
#define UTILS_H

#include "BasicDataStructure.h"
#include <vector>
#include <iostream>

using namespace std;


// ! Matrix operation
double MatrixDet(double dxdt[2][2]);

double MatrixDet(double dxdt[3][3]);

void Matrix2DInverse(double dxdt[2][2], double dtdx[2][2]);

void Matrix3DInverse(double dxdt[3][3], double dtdx[3][3]);

// ! Debug Function
void DebugMatVector(vector<vector<double>> mat_debug, int e, string fname);

void DebugVisualizeMat(Mat mat_debug, string fname);

void DebugVisualizeVec(Vec vec_debug, string fname);

void DebugVisualizeIS(IS is_debug, string fname);

void DebugVisualizeMatSelf(Mat mat_debug, string fname);

void DebugVisualizeVecSelf(Vec vec_debug, string fname);

void DebugVisualizeISSelf(IS is_debug, string fname);

#endif
