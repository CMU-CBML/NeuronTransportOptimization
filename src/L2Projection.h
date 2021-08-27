#ifndef L2PROJECTION_H
#define L2PROJECTION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
#include "Utils.h"
#include <cmath>
using namespace std;

//problem setting
class L2Projection
{
private:
    static const bool ReadIC = false;
    static const bool ReadBC = false;
    // static const bool ReadDesire = true;
    static const bool ReadDesire = false;

    static const bool VisualizeIC = true;
    static const bool VisualizeBC = true;
    static const bool VisualizeDesire = true;

    static const int phy_dim = 3;
    static const int dim = 2;
    static const int degree = 3;
    static const int bzpt_num = 16;
    static const int bzpt_num_3d = 64;

    // * Neuron Model equation
    // const int state_num = 7;
    // const int ctrl_num = 4;
    // const int result_num = 7;
    // * Burger's equation
    // const int state_num = 2;
    // const int ctrl_num = 2;
    // const int result_num = 6;
    // * Diffusion equation and Convection-diffusion equation
    static const int state_num = 1;
    static const int ctrl_num = 1;
    static const int result_num = 3;

    static const int time_int = 0; // * 0 - steady state; 1 - trapezoidal; 2 - rectangle
public:
    int n_bzmesh;
    int n_bcpt;
    int n_bcval;

    vector<vector<int>> ele_process;
    vector<Element2D> bzmesh_process;
    vector<Element3D> bzmesh_process_3d;

    DM dm_bzmesh;

    vector<double> var;
    vector<Vertex2D> pts;
    vector<Element2D> mesh;

    vector<Vertex3D> pts_3d;
    vector<Element3D> mesh_3d;

    vector<double> val_ini[2];
    vector<double> val_bc[2];
    vector<double> val_desire[2];
    vector<int> bc_flag;       // global node index -> index without bc pts
    vector<int> nonbc_mapping; // index without bc pts -> global node index
    vector<int> bc_mapping;    // index of bc pts -> global node index

    string output_dir;
    string work_dir;

private:
    PetscErrorCode ierr;
    MPI_Comm comm;
    int mpiErr;
    int comRank;
    int comSize;
    int nProcess;
    KSP ksp_L2;
    PC pc_L2;
    Mat M;
    Vec X, Rhs;
    int nTstep, nPoint;

    vector<double> Gpt;
    vector<double> wght;
    const double PI = 4 * atan(1.0);

    void GaussInfo(int ng);
    void BasisFunction2D(double u, double v, int nen, const vector<array<double, phy_dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, double &detJ);
    void BasisFunction3D(double u, double v, double w, int nen, const vector<array<double, phy_dim>> &pt, const vector<array<double, bzpt_num_3d>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, double &detJ);
    void VecAssembly(vector<double> Eqtmp, const vector<int> &IEN, Vec &Gqtmp);
    void MatrixAssembly(vector<vector<double>> Emat, const vector<int> &IEN, Mat &Gmat);
    void ComputePhysCoor(const vector<Vertex2D> &cpts, const vector<int> &IEN, vector<double> &Nx, vector<double> &PhysCoor);
    void ComputePhysCoor(const vector<Vertex3D> &cpts, const vector<int> &IEN, vector<double> &Nx, vector<double> &PhysCoor);
    void ComputeRhsVec(double &RhsVal, vector<double> &Nx, const double detJ, vector<double> &RhsVec);
    void ComputeMassMatrix(vector<double> &Nx, const double detJ, vector<vector<double>> &MassMat);
    void ComputeDesireStateFunction(const vector<double> &PhysCoor, double t, double &result);

    void BuildLinearSystemProcess2D(const vector<Vertex2D> &cpts);
    void BuildLinearSystemProcess3D(const vector<Vertex3D> &cpts);
    void CollectAndUpdateResult(Vec dx);

private:
    void SetVariables(string fn_par);
    void SetDesireState(string fn_in, string fn_out);
    void SetInitialCondition(string fn_in, string fn_out);
    void SetBoundaryCondition(string fn_in, string fn_out);
    void SetBoundaryMapping();

    void TXTWriteIC(string fn_out);
    void TXTWriteBC(string fn_out);
    void TXTWriteDesire(string fn_out);

    void VTKVisualizeIC(string fn_out);
    void VTKVisualizeBC(string fn_out);
    void VTKVisualizeDesire(string fn_out);

    void ReadMesh(string fn);
    void CreateDMPlexMesh(string fn);
    void ReadBezierElementProcess(string fn);
    void ReadVelocityField(string fn);
    void AssignProcessor(string fn);

public:
    L2Projection();
    ~L2Projection();
    bool GetReadDesire() const;

    void InitializeUserSetting(string fn_in, string fn_out);

    void DesireStateFunction(double x, double y, double z, double t, double result[state_num]) const;
    void RunL2Projection(string fn);
};
#endif