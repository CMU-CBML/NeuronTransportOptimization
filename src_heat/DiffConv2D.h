#ifndef DiffConv2D_H
#define DiffConv2D_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "BasicDataStructure.h"
#include "UserSetting2D.h"
#include "time.h"

using namespace std;



class DiffConv2D
{
private:
    vector<double> Gpt;
    vector<double> wght;
    const double PI = 4 * atan(1.0);

    PetscErrorCode ierr;
    MPI_Comm comm;
    int mpiErr;
    int comRank;
    int comSize;
    int nProcess;
    string work_dir;

    int rstart, rend;
    int n_bzmesh;
    vector<int> ele_process;
    vector<Element2D> bzmesh_process;

    SNES snes;
    KSP ksp;
    PC pc;

    double dt;
    double velocity_max;
    double nu;
    double rou;
    double alphaM;
    double alphaF;
    double Gama;
    double fx, fy, fz; //elemental body force

    // Scale of the problem
    int nTstep, nPoint;
    int n_bcval; // Number of boundary vals we set, depends on the problem

    // constant parameters
    vector<double> par; //Dn0, v_plus, v_minus, k+, k-,k'+,k'-,c0,mu
    double alpha0, alpha1, alpha2;
    double beta1, beta2;

    // state variables
    vector<double> State[state_num];
    // control variables
    vector<double> Ctrl[ctrl_num];
    // penalty variables
    vector<double> Lambda[state_num];

    // Unit Np * Np Matrix
    Mat M, K, P[dim], Conv, Stable;

    // KKT system matrix and vector
    Mat Asubmat[9];
    Mat PCsubmat[9]; // For Preconditioner
    Vec bsub[3];
    IS isg[3];

    // Variables for solving iteration
    Vec Y_d, Y_bc;
    Vec Y_ini, U_ini, L_ini; // Initial value in time space
    Vec Y_k, U_k, L_k;       // Initial value in iteration
    Vec Res_Desire, Res_nl;

    Vec dX;
    Vec X, ResVec;
    Mat TanMat_tmp, TanMat;
    Mat PCMat_tmp, PCMat;

    vector<double> Vel;
    vector<double> Pre;
    vector<double> gh; //prescribed boundary velocities

public:
    DiffConv2D();
    ~DiffConv2D();

private:
    void BodyForce(double x, double y, double &Fx, double &Fy);

    void Debug();

    /*MPI computing*/
    void ReadBezierElementProcess(string fn);
    void AssignProcessor(const UserSetting2D *ctx);

    /*Analysis*/
    void GaussInfo(int ng);
    void InitializeProblem(const UserSetting2D *ctx);
    void BasisFunction(double u, double v, int nen, const vector<array<double, dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, double dudx[dim][dim], double &detJ);
    void BasisFunctionCoarseSpace(double u, double v, int nen, const vector<array<double, dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, double &detJ);

    // void BasisFunction(double u, double v, double w, int nen, const vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat, vector<double> &Nx, vector<array<double, 3>> &dNdx, vector<array<array<double, 3>, 3>> &dN2dx2, double dudx[3][3], double& detJ);
    void PointFormValue(vector<double> &Nx, const vector<double> &U, double Value);
    void PointFormGrad(vector<array<double, dim>> &dNdx, const vector<double> &U, double Value[dim]);
    void PointFormHess(vector<array<array<double, dim>, dim>> &d2Ndx2, const vector<double> &U, double Value[dim][dim]);

    void ComputeMassMatrix(vector<double> &Nx, const double detJ, vector<vector<double>> &MassMat);
    void ComputeStiffMatrix(vector<array<double, dim>> &dNdx, const double detJ, vector<vector<double>> &StiffMat);
    void ComputeConvectMatrix(vector<double> &v, vector<double> &Nx, vector<array<double, dim>> &dNdx, const double detJ, vector<vector<double>> &ConvectMat);
    void ComputeParMatrix(vector<double> &Nx, vector<array<double, dim>> &dNdx, const double detJ, int dir, vector<vector<double>> &ParMat);

    void LocalL2Projection(int e, vector<double> w, vector<vector<double>> &u_proj);

    void ComputeStableMatrix(const vector<double> val_ini[state_num], double h, vector<double> wdNdx_proj, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, const double detJ, const vector<int> &IEN, vector<vector<double>> &StableMat);
    void ComputeStableMatrix_Convection(const vector<double> w, double h, vector<double> wdNdx_proj, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<array<array<double, dim>, dim>> &dN2dx2, const double detJ, const vector<int> &IEN, vector<vector<double>> &StableMat);
    void ComputeStableMatrix_Convection2(const vector<double> w, double h, vector<vector<double>> wdNdx_proj, vector<double> &Nx, vector<array<double, dim>> &dNdx, vector<double> &Nx_disc, vector<array<double, dim>> &dNdx_disc, const double detJ, const double detJ_disc, const vector<int> &IEN, vector<vector<double>> &StableMat);

    //void ComputeConvectionMatrix(vector<double>& Nx, vector<array<double, 3>>& dNdx, const double detJ, const vector<double> &U, vector<vector<double>> ConvectMat);
    void ComputeResVector(const vector<int> bc_flag, const vector<double> w, const vector<double> val_ini[state_num], vector<double> &Nx, vector<array<double, dim>> &dNdx, const vector<int> &IEN, const double detJ, vector<double> &qtmp);

    void ResnlAssembly(vector<double> Eqtmp, const vector<int> &IEN, Vec &Gqtmp);
    void StableMatAssembly(vector<vector<double>> Emat, const vector<int> &IEN, Mat &Gmat);
    void MatrixAssembly(vector<vector<double>> Emat, const vector<int> &IEN, Mat &Gmat);
    void BuildLinearSystemProcess(const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num], const vector<double> val_ini[state_num]);
    void BuildLinearSystemProcessAsubmat(const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num], const vector<double> val_ini[state_num]);

    void BuildResVectorProcess(const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num], const vector<double> val_ini[state_num], const vector<int> bc_flag);

    void GetMatrixPosition(int row, int n_var, int &i_point, int &i_var, int &i_tstep);

    void ElementFormMatrixA11(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA12(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA13(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA21(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA22(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA23(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA31(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA32(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);
    void ElementFormMatrixA33(const vector<int> &IEN, vector<vector<double>> EM, vector<vector<double>> EK, vector<vector<double>> EConv, vector<vector<double>> EPx, vector<vector<double>> EPy, vector<vector<double>> EStable);

    void FormMatrixA11(Mat M, Mat K, Mat &A);
    void FormMatrixA12(Mat M, Mat K, Mat &A);
    void FormMatrixA13(Mat M, Mat K, Mat P[dim], Mat &A);
    void FormMatrixA21(Mat M, Mat K, Mat &A);
    void FormMatrixA22(Mat M, Mat K, Mat &A);
    void FormMatrixA23(Mat M, Mat K, Mat &A);
    void FormMatrixA31(Mat M, Mat K, Mat P[dim], Mat &A);
    void FormMatrixA32(Mat M, Mat K, Mat &A);
    void FormMatrixA33(Mat M, Mat K, Mat &A);

    void AssembleGlobalMatrix(Mat &A);
    void TangentMatSetup();
    void TangentMatSetupAsubmat();
    void PCMatSetup();

    void FormDesireVec(const UserSetting2D *ctx);
    void FormResVecb1(Vec x);
    void FormResVecb2(Vec x);
    void FormResVecb3(Vec x);
    void ResidualVecSetup(Vec x);

    void ApplyBoundaryCondition(const UserSetting2D *ctx, int flag);

    void CollectAndUpdateResult(Vec dx, Vec x);

    PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx);
    PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx);

    void DebugSubMat();

    /*Postprocessing*/
    void ResultCal_Bezier(double u, double v, int time, const Element2D &bzel, double pt[2], double result[result_num], double dudx[3], double &detJ);
    void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 4>> sele, int time, int step, string fn);
    void VisualizeVTK_PhysicalDomain(int time, int step, string fn);
    void VisualizeVTK_ControlMesh(const vector<Vertex2D> &pts, const vector<Element2D> &mesh, int step, string fn);
    void VisualizeVTK_ControlMesh_Burger(const vector<Vertex2D> &pts, const vector<Element2D> &mesh, int time, int step, string fn);
    void VisualizeVTK_ControlMesh_Heat(const vector<Vertex2D> &pts, const vector<Element2D> &mesh, int time, int step, string fn);

    //void ResultCal_Bezier(double u, double v, const Element2D& bzel, double pt[3], double result[4], double dudx[3], double& detJ);

public:
    /*Preprocessing*/
    void InitializeProblem(const int ndof, const int n_bz, const vector<double> &Vel0, const vector<double> &Pre0, const vector<double> &var);
    void AssignProcessor(vector<vector<int>> &ele_proc);
    void Run(const vector<Vertex2D> &cpts, const vector<Element2D> &tmesh, const vector<array<double, 2>> &velocity_bc, string fn);

    void Run(const UserSetting2D *ctx);
    void Run_SNES(const UserSetting2D *ctx);
};

#endif