#include "L2Projection.h"

L2Projection::L2Projection()
{
    comm = MPI_COMM_WORLD;
    mpiErr = MPI_Comm_rank(comm, &comRank);
    mpiErr = MPI_Comm_size(comm, &comSize);
    nProcess = comSize;
}

L2Projection::~L2Projection()
{
}

void L2Projection::AssignProcessor(vector<vector<int>> &ele_proc)
{
    /*Assign the partitioned bezier elements to this processor */
    for (int i = 0; i < ele_proc[comRank].size(); i++)
        ele_process.push_back(ele_proc[comRank][i]);
}

void L2Projection::ReadBezierElementProcess(string fn)
{
    string stmp;
    int itmp, itmp1;
    int npts, neles, nfunctions;

    string fname_cmat = fn + "cmat.txt";
    ifstream fin_cmat;
    fin_cmat.open(fname_cmat);
    if (fin_cmat.is_open())
    {
        fin_cmat >> neles;
        // cout << neles << endl;
        bzmesh_process.resize(ele_process.size());

        int add(0);
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

    MPI_Barrier(comm);

    string fname_bzpt = fn + "bzpt.txt";
    ifstream fin_bzpt;
    fin_bzpt.open(fname_bzpt);
    // cout << fname_bzpt <<endl;

    if (fin_bzpt.is_open())
    {
        int add = 0;
        fin_bzpt >> npts;
        // cout << npts << endl;
        getline(fin_bzpt, stmp);
        for (int e = 0; e < neles; e++)
        {
            // if (e == ele_process[add])
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
    MPI_Barrier(comm);
}

void L2Projection::GaussInfo(int ng)
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

void L2Projection::BasisFunction2D(double u, double v, int nen, const vector<array<double, phy_dim>> &pt, const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx, vector<array<double, dim>> &dNdx, double &detJ)
{
    double Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
    double Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};

    double dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
    double dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};

    double dNdt[bzpt_num][dim];
    double dN2dt2[bzpt_num][dim][dim];
    double Nx_bz[bzpt_num];
    double dNdx_bz[bzpt_num][dim];

    Nx.clear();
    dNdx.clear();
    Nx.resize(nen, 0);
    dNdx.resize(nen, {0});

    int i, j, k, a, b, c, loc;
    loc = 0;

    for (j = 0; j < 4; j++)
    {
        for (k = 0; k < 4; k++)
        {
            Nx_bz[loc] = Nu[k] * Nv[j];
            dNdt[loc][0] = dNdu[k] * Nv[j];
            dNdt[loc][1] = Nu[k] * dNdv[j];
            loc++;
        }
    }

    double dxdt[2][2] = {{0}};
    for (loc = 0; loc < bzpt_num; loc++)
        for (a = 0; a < 2; a++)
            for (b = 0; b < 2; b++)
                dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

    double dtdx[2][2] = {{0}};
    Matrix2DInverse(dxdt, dtdx);

    //1st derivatives
    for (i = 0; i < bzpt_num; i++)
    {
        dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
        dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
    }

    detJ = MatrixDet(dxdt);
    detJ = 0.25 * detJ;

    for (i = 0; i < nen; i++)
    {
        for (j = 0; j < bzpt_num; j++)
        {
            Nx[i] += cmat[i][j] * Nx_bz[j];
            for (int m = 0; m < 2; m++)
            {
                dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
            }
        }
    }
}

void L2Projection::BasisFunction3D(double u, double v, double w, int nen,
                                   const vector<array<double, phy_dim>> &pt,
                                   const vector<array<double, bzpt_num_3d>> &cmat,
                                   vector<double> &Nx,
                                   vector<array<double, dim>> &dNdx, double &detJ)
{
    double Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u, 3. * (1. - u) * u * u, u * u * u};
    double Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v, 3. * (1. - v) * v * v, v * v * v};
    double Nw[4] = {(1. - w) * (1. - w) * (1. - w), 3. * (1. - w) * (1. - w) * w, 3. * (1. - w) * w * w, w * w * w};

    double dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u, 3. * (2. - 3. * u) * u, 3. * u * u};
    double dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v, 3. * (2. - 3. * v) * v, 3. * v * v};
    double dNdw[4] = {-3. * (1. - w) * (1. - w), 3. - 12. * w + 9. * w * w, 3. * (2. - 3. * w) * w, 3. * w * w};

    double dN2du2[4] = {6. * (1. - u), -12. + 18. * u, 6. - 18. * u, 6. * u};
    double dN2dv2[4] = {6. * (1. - v), -12. + 18. * v, 6. - 18. * v, 6. * v};
    double dN2dw2[4] = {6. * (1. - w), -12. + 18. * w, 6. - 18. * w, 6. * w};

    double dNdt[bzpt_num][3];
    double dN2dt2[bzpt_num][3][3];
    double Nx_bz[bzpt_num];
    double dNdx_bz[bzpt_num][3];
    double dN2dx2_bz[bzpt_num][3][3];

    Nx.clear();
    dNdx.clear();
    Nx.resize(nen, 0);
    dNdx.resize(nen, {0});

    int i, j, k, a, b, c, loc;
    loc = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            for (k = 0; k < 4; k++)
            {
                Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
                dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
                dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
                dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
                dN2dt2[loc][0][0] = dN2du2[k] * Nv[j] * Nw[i];
                dN2dt2[loc][0][1] = dNdu[k] * dNdv[j] * Nw[i];
                dN2dt2[loc][0][2] = dNdu[k] * Nv[j] * dNdw[i];
                dN2dt2[loc][1][0] = dNdu[k] * dNdv[j] * Nw[i];
                dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j] * Nw[i];
                dN2dt2[loc][1][2] = Nu[k] * dNdv[j] * dNdw[i];
                dN2dt2[loc][2][0] = dNdu[k] * Nv[j] * dNdw[i];
                dN2dt2[loc][2][1] = Nu[k] * dNdv[j] * dNdw[i];
                dN2dt2[loc][2][2] = Nu[k] * Nv[j] * dN2dw2[i];
                loc++;
            }
        }
    }

    double dxdt[3][3] = {{0}};
    for (loc = 0; loc < bzpt_num; loc++)
        for (a = 0; a < 3; a++)
            for (b = 0; b < 3; b++)
                dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

    double dtdx[3][3] = {{0}};
    Matrix3DInverse(dxdt, dtdx);

    //1st derivatives
    for (i = 0; i < bzpt_num; i++)
    {
        dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0] + dNdt[i][2] * dtdx[2][0];
        dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1] + dNdt[i][2] * dtdx[2][1];
        dNdx_bz[i][2] = dNdt[i][0] * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
    }

    detJ = MatrixDet(dxdt);
    detJ = 0.125 * detJ;

    //2nd derivatives
    double dx2dt2[3][9] = {{0}};
    double dt2dx2[3][9] = {{0}};
    for (int l = 0; l < 3; l++)
        for (loc = 0; loc < bzpt_num; loc++)
            for (a = 0; a < 3; a++)
                for (b = 0; b < 3; b++)
                    dx2dt2[l][3 * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                for (a = 0; a < 3; a++)
                    for (b = 0; b < 3; b++)
                        for (c = 0; c < 3; c++)
                            dt2dx2[c][3 * i + j] -= dx2dt2[k][3 * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

    for (loc = 0; loc < bzpt_num; loc++)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                dN2dx2_bz[loc][i][j] = 0.;

    for (loc = 0; loc < bzpt_num; loc++)
    {
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                for (a = 0; a < 3; a++)
                {
                    for (b = 0; b < 3; b++)
                    {
                        dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
                    }
                    dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][3 * i + j];
                }
            }
        }
    }
}

void L2Projection::VecAssembly(vector<double> Eqtmp, const vector<int> &IEN, Vec &Gqtmp)
{
    int a, b, nen = IEN.size();
    PetscInt *rows;
    PetscReal *vals;

    PetscMalloc1(IEN.size() * state_num * nTstep, &rows);
    PetscMalloc1(IEN.size() * state_num * nTstep, &vals);

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < nen; i++)
        {
            for (int k = 0; k < state_num; k++)
            {
                int A = IEN[i] + k * nPoint + j * nPoint * state_num;
                rows[i + nen * k + j * nen * state_num] = A;
                vals[i + nen * k + j * nen * state_num] = Eqtmp[i + nen * k + j * nen * state_num];
            }
        }
    }
    VecSetValues(Gqtmp, nen * state_num * nTstep, rows, vals, ADD_VALUES);

    PetscFree(rows);
    PetscFree(vals);

    delete rows, vals;
}

void L2Projection::MatrixAssembly(vector<vector<double>> Emat, const vector<int> &IEN, Mat &Gmat)
{
    int i, j, A, B, m, n;
    int row_start, row_end, row_now;
    int add = 0;

    PetscInt *nodeList = new PetscInt[IEN.size()];
    PetscReal *tmpGmat = new PetscReal[IEN.size() * IEN.size()];

    for (m = 0; m < IEN.size(); m++)
    {
        A = IEN[m];
        nodeList[m] = A;
        for (n = 0; n < IEN.size(); n++)
        {
            B = IEN[n];
            tmpGmat[add] = Emat[m][n];
            add++;
        }
    }

    MatSetValues(Gmat, IEN.size(), nodeList, IEN.size(), nodeList, tmpGmat, ADD_VALUES);
    delete nodeList;
    delete tmpGmat;
}

void L2Projection::ComputePhysCoor(const vector<Vertex2D> &cpts, const vector<int> &IEN, vector<double> &Nx, vector<double> &PhysCoor)
{
    PhysCoor.clear();
    PhysCoor.resize(dim, 0);
    int i, j;
    for (i = 0; i < IEN.size(); i++)
    {
        for (j = 0; j < dim; j++)
        {
            PhysCoor[j] += Nx[i] * cpts[IEN[i]].coor[j];
        }
    }
}

void L2Projection::ComputePhysCoor(const vector<Vertex3D> &cpts, const vector<int> &IEN, vector<double> &Nx, vector<double> &PhysCoor)
{
    PhysCoor.clear();
    PhysCoor.resize(dim, 0);
    int i, j;
    for (i = 0; i < IEN.size(); i++)
    {
        for (j = 0; j < dim; j++)
        {
            PhysCoor[j] += Nx[i] * cpts[IEN[i]].coor[j];
        }
    }
}

void L2Projection::ComputeDesireStateFunction(const vector<double> &PhysCoor, double t, double &result)
{
    double eps = 1e-6;
    double x0 = PhysCoor[0];
    double x1 = PhysCoor[1];

    if (abs(x0) <= eps)
    {
        // result = pow((2.0 * x0 - 1.0), 2.0) * pow((2.0 * x1 - 1.0), 2.0); // test1 and test2
        result = -160.0 * (x1 * x1) + 40.0; // test3 and test4
        // cout << "x0: " << x0 << " x1: " << x1 << " result: " << result << endl;
    }
    else
    {
        result = 0.0;
    }
}

void L2Projection::ComputeRhsVec(double &RhsVal, vector<double> &Nx, const double detJ, vector<double> &RhsVec)
{
    int nen = Nx.size();
    for (int i = 0; i < nen; i++)
        RhsVec[i] += RhsVal * Nx[i] * detJ;
}

void L2Projection::ComputeMassMatrix(vector<double> &Nx, const double detJ, vector<vector<double>> &MassMat)
{
    int a, b, nen = Nx.size();
    for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
            MassMat[a][b] += Nx[a] * Nx[b] * detJ;
}

void L2Projection::BuildLinearSystemProcess2D(const vector<Vertex2D> &cpts)
{
    //Build linear system in each process
    int e;
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for " << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        double h;

        int nen, A;
        double ux_bc, uy_bc, uz_bc, p_bc;

        double detJ;
        vector<double> Nx;
        vector<array<double, dim>> dNdx;
        vector<vector<double>> Mtmp;
        vector<double> Rhstmp;

        vector<double> PhysCoor;
        double DesireVal = 0;
        double BCVal = 0;
        double ICVal = 0;

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0);
        dNdx.clear();
        dNdx.resize(nen, {0});

        Rhstmp.clear();
        Rhstmp.resize(nen, 0);

        Mtmp.clear();
        Mtmp.resize(nen);
        for (int i = 0; i < nen; i++)
        {
            Mtmp[i].resize(nen, 0.);
        }
        // getchar();

        int debug_count = 0;
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                BasisFunction2D(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                ComputeMassMatrix(Nx, detJ, Mtmp);
                ComputePhysCoor(cpts, bzmesh_process[e].IEN, Nx, PhysCoor);
                ComputeDesireStateFunction(PhysCoor, 0, DesireVal);

                ComputeRhsVec(DesireVal, Nx, detJ, Rhstmp);
                //ComputeRhsVec(BCVal, Nx, detJ, Rhstmp_bc);
                // ComputeRhsVec(ICVal, Nx, detJ, Rhstmp_ini);
                debug_count++;
            }
        }
        // MPI_Barrier(comm);

        //Start element matrix assembly
        VecAssembly(Rhstmp, bzmesh_process[e].IEN, Rhs);
        // VecAssembly(Rhstmp_ini, bzmesh_process[e].IEN, Rhs_ini);
        MatrixAssembly(Mtmp, bzmesh_process[e].IEN, M);
    }
    cout << "Process " << comRank << " :complete build matrix and vector!\n";
    VecAssemblyBegin(Rhs);
    // VecAssemblyBegin(Rhs_ini);
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);

    MPI_Barrier(comm);
    // VecAssemblyBegin(Res_nl);
}

void L2Projection::BuildLinearSystemProcess3D(const vector<Vertex3D> &cpts)
{
    //Build linear system in each process
    int e;
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for " << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        double h;

        int nen, A;
        double ux_bc, uy_bc, uz_bc, p_bc;

        double detJ;
        vector<double> Nx;
        vector<array<double, dim>> dNdx;
        vector<vector<double>> Mtmp;
        vector<double> Rhstmp;

        vector<double> PhysCoor;
        double DesireVal = 0;
        double BCVal = 0;
        double ICVal = 0;

        nen = bzmesh_process_3d[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0);
        dNdx.clear();
        dNdx.resize(nen, {0});

        Rhstmp.clear();
        Rhstmp.resize(nen, 0);

        Mtmp.clear();
        Mtmp.resize(nen);
        for (int i = 0; i < nen; i++)
        {
            Mtmp[i].resize(nen, 0.);
        }
        // getchar();

        int debug_count = 0;
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                for (int k = 0; k < Gpt.size(); k++)
                {
                    BasisFunction3D(Gpt[i], Gpt[j], Gpt[k], nen, bzmesh_process_3d[e].pts, bzmesh_process_3d[e].cmat, Nx, dNdx, detJ);
                    detJ = wght[i] * wght[j] * detJ;
                    ComputeMassMatrix(Nx, detJ, Mtmp);
                    ComputePhysCoor(cpts, bzmesh_process_3d[e].IEN, Nx, PhysCoor);
                    ComputeDesireStateFunction(PhysCoor, 0, DesireVal);

                    ComputeRhsVec(DesireVal, Nx, detJ, Rhstmp);
                    //ComputeRhsVec(BCVal, Nx, detJ, Rhstmp_bc);
                    // ComputeRhsVec(ICVal, Nx, detJ, Rhstmp_ini);
                    debug_count++;
                }
            }
        }
        // MPI_Barrier(comm);

        //DebugMatVector(Stabletmp, e, work_dir + "debug/Stabletmp.txt");

        //Start element matrix assembly
        VecAssembly(Rhstmp, bzmesh_process_3d[e].IEN, Rhs);
        // VecAssembly(Rhstmp_ini, bzmesh_process[e].IEN, Rhs_ini);
        MatrixAssembly(Mtmp, bzmesh_process_3d[e].IEN, M);
    }
    cout << "Process " << comRank << " :complete build matrix and vector!\n";
    VecAssemblyBegin(Rhs);
    // VecAssemblyBegin(Rhs_ini);
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);

    MPI_Barrier(comm);
    // VecAssemblyBegin(Res_nl);
}

void L2Projection::CollectAndUpdateResult(Vec dx)
{
    Vec dx_seq;
    VecScatter scatter_ctx;
    PetscReal *dx_array;
    VecScatterCreateToAll(dx, &scatter_ctx, &dx_seq);
    VecScatterBegin(scatter_ctx, dx, dx_seq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter_ctx, dx, dx_seq, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(dx_seq, &dx_array);
    MPI_Barrier(comm);

    // ! for linear system

    for (uint i = 0; i < 2; i++)
        val_desire[i].resize(pts.size());

    for (uint j = 0; j < nTstep; j++)
    {
        for (uint i = 0; i < nPoint; i++)
        {
            int A, B;
            for (uint k = 0; k < state_num; k++)
            {
                A = i + j * nPoint;
                B = i + k * nPoint + j * nPoint * state_num;
                val_desire[k][A] = PetscRealPart(dx_array[B]);
            }
        }
    }

    VecRestoreArray(dx_seq, &dx_array);
    VecScatterDestroy(&scatter_ctx);
    VecDestroy(&dx_seq);
}

void L2Projection::RunL2Projection(const vector<Vertex2D> &cpts, const vector<Element2D> &tmesh, string fn)
{

    PetscPrintf(PETSC_COMM_WORLD, "Start L2 Projection...\n");
    work_dir = fn;
    //! L2-projection of desire function

    ReadBezierElementProcess(work_dir);
    pts = cpts;
    mesh = tmesh;

    nPoint = pts.size();
    PetscPrintf(PETSC_COMM_WORLD, "nPoints:%d\n", nPoint);

    ierr = MatCreate(PETSC_COMM_WORLD, &M);
    ierr = MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(M, MATMPIAIJ);
    if (dim == 2)
        ierr = MatMPIAIJSetPreallocation(M, 64, NULL, 64, NULL);
    else
        ierr = MatMPIAIJSetPreallocation(M, 125, NULL, 125, NULL);
    ierr = MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(M);
    ierr = PetscObjectSetName((PetscObject)M, "M");

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &X);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, nPoint * state_num * nTstep, &Rhs);

    ierr = VecSet(Rhs, 0.0);

    GaussInfo(4);
    PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");

    time_t t0, t1;

    t0 = time(NULL);
    if (dim == 2)
        BuildLinearSystemProcess2D(pts);
    else
        BuildLinearSystemProcess3D(pts_3d);
    t1 = time(NULL);

    VecAssemblyEnd(Rhs);
    PetscPrintf(PETSC_COMM_WORLD, "Done Vector Assembly...\n");
    DebugVisualizeVec(Rhs, work_dir + "debug/Rhs_L2.m");

    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    DebugVisualizeMat(M, work_dir + "debug/M_L2.m");
    PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 - t0);
    // //VecAssemblyEnd(Rhs_ini);
    KSPCreate(PETSC_COMM_WORLD, &ksp_L2);
    KSPSetOperators(ksp_L2, M, M);

    KSPSolve(ksp_L2, Rhs, X);
    DebugVisualizeVec(X, work_dir + "debug/X_L2.m");
    CollectAndUpdateResult(X);

    // VisualizeVTK_ControlMesh_Heat(pts, mesh, 0, 0, work_dir);
    // VisualizeVTK_PhysicalDomain(0, 0, work_dir);
    //KSPSolve(ksp, Rhs_ini, Y_ini);
}

void L2Projection::VisualizeVTK_ControlMesh_Heat(const vector<Vertex2D> &spt,
                                                 const vector<Element2D> &mesh,
                                                 int time, int step, string fn)
{
    stringstream ss;
    ss << step << "_" << time;

    string fname = fn + "L2Result_CM_" + ss.str() + ".vtk";
    // string fname = fn + "controlmesh_VelocityPressure.vtk";
    ofstream fout;
    fout.open(fname.c_str());
    unsigned int i;
    if (fout.is_open())
    {
        fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET "
                "UNSTRUCTURED_GRID\n";
        fout << "POINTS " << spt.size() << " float\n";
        for (i = 0; i < spt.size(); i++)
        {
            fout << std::setprecision(9) << spt[i].coor[0] << std::fixed << " "
                 << std::setprecision(9) << spt[i].coor[1] << std::fixed << " "
                 << std::setprecision(9) << spt[i].coor[2] << std::fixed << "\n";
        }
        fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
        for (i = 0; i < mesh.size(); i++)
        {
            fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " "
                 << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
        }
        fout << "\nCELL_TYPES " << mesh.size() << "\n";
        for (i = 0; i < mesh.size(); i++)
        {
            fout << "9\n";
        }
        fout << "POINT_DATA " << spt.size() << "\n";
        fout << "VECTORS y float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << 0 << " " << std::setprecision(9) << 0
                 << " " << val_desire[0][i + time * nPoint] << "\n";
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }
}

void L2Projection::VisualizeVTK_PhysicalDomain(int time, int step, string fn)
{
    vector<array<double, 3>> spt_all; // sample points
    vector<double> sresult_all;
    vector<array<int, 4>> sele_all;

    const int ns_ele = 9;
    const int ns_pt = 16;
    double detJ;
    int num_bzmesh_ele = bzmesh_process.size();
    double spt_proc[num_bzmesh_ele * ns_pt * 3];
    double sresult_proc[num_bzmesh_ele * ns_pt * result_num];
    int sele_proc[num_bzmesh_ele * ns_ele * 4];
    for (unsigned int e = 0; e < num_bzmesh_ele; e++)
    {
        int ns(4);
        vector<double> su;
        su.clear();
        su.resize(ns);
        for (int i = 0; i < ns; i++)
        {
            su[i] = double(i) / (double(ns) - 1.);
        }

        int loc(0);

        // * Compute cordinates and nodal results on sample points
        for (int a = 0; a < ns; a++)
        {
            for (int b = 0; b < ns; b++)
            {
                double pt1[3], dudx[3];
                double result[result_num];
                ResultCal_Bezier(su[b], su[a], time, bzmesh_process[e], pt1, result,
                                 dudx, detJ);
                // cout << " spt: " << loc << " " << pt1[0] << " " << pt1[1] << " " <<
                // pt1[2] << "\n";
                spt_proc[e * ns_pt * 3 + loc * 3 + 0] = pt1[0];
                spt_proc[e * ns_pt * 3 + loc * 3 + 1] = pt1[1];
                spt_proc[e * ns_pt * 3 + loc * 3 + 2] = pt1[2]; // ! for 2D problem
                for (int c = 0; c < result_num; c++)
                {
                    sresult_proc[result_num * ns_pt * e + loc * result_num + c] =
                        result[c];
                }
                loc++;
            }
        }

        // getchar();
        // * Compute mesh connectivity of sample elements
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
    MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0,
               PETSC_COMM_WORLD);
    MPI_Barrier(comm);

    if (comRank == 0)
    {
        spts = (double *)malloc(sizeof(double) * 3 * n_bzmesh * ns_pt);
        sresults = (double *)malloc(sizeof(double) * result_num * n_bzmesh * ns_pt);
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
                displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * result_num * ns_pt;
            displs_seles[i] =
                displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4 * ns_ele;
        }

        for (int i = 0; i < nProcess; i++)
        {
            recvcounts_spts[i] = num_bzmesh_eles[i] * ns_pt * 3;
            recvcounts_sresults[i] = num_bzmesh_eles[i] * ns_pt * result_num;
            recvcounts_seles[i] = num_bzmesh_eles[i] * ns_ele * 4;
        }
    }

    MPI_Gatherv(spt_proc, num_bzmesh_ele * ns_pt * 3, MPI_DOUBLE, spts,
                recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Gatherv(sresult_proc, num_bzmesh_ele * ns_pt * result_num, MPI_DOUBLE,
                sresults, recvcounts_sresults, displs_sresults, MPI_DOUBLE, 0,
                PETSC_COMM_WORLD);
    MPI_Gatherv(sele_proc, num_bzmesh_ele * ns_ele * 4, MPI_INT, seles,
                recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

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
                for (int c = 0; c < result_num; c++)
                {
                    sresult_all.push_back(
                        sresults[i * result_num * ns_pt + j * result_num + c]);
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
        WriteVTK(spt_all, sresult_all, sele_all, time, step, fn);
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

void L2Projection::WriteVTK(const vector<array<double, 3>> spt,
                            const vector<double> sdisp,
                            const vector<array<int, 4>> sele, int time, int step,
                            string fn)
{
    stringstream ss;
    ss << step << "_" << time;

    string fname = fn + +"L2Result_physical_" + ss.str() + ".vtk";
    ofstream fout;
    fout.open(fname.c_str());
    unsigned int i;
    if (fout.is_open())
    {
        fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET "
                "UNSTRUCTURED_GRID\n";
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
            fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2]
                 << " " << sele[i][3] << '\n';
        }
        fout << "\nCELL_TYPES " << sele.size() << '\n';
        for (i = 0; i < sele.size(); i++)
        {
            fout << "9\n";
        }
        fout << "POINT_DATA " << sdisp.size() / result_num << "\n";
        fout << "VECTORS y float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << sdisp[i * result_num + 0] << std::fixed
                 << " " << 0 << std::fixed << " " << 0 << "\n";
        fout << "VECTORS u float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << sdisp[i * result_num + 1] << std::fixed
                 << " " << 0 << std::fixed << " " << 0 << "\n";
        fout << "VECTORS lambda float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << sdisp[i * result_num + 2] << std::fixed
                 << " " << 0 << std::fixed << " " << 0 << "\n";
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }
}

void L2Projection::ResultCal_Bezier(double u, double v, int time,
                                    const Element2D &bzel, double pt[2],
                                    double result[result_num], double dudx[3],
                                    double &detJ)
{
    double dUdx[dim][dim];
    vector<double> Nx(bzel.IEN.size());
    vector<array<double, dim>> dNdx(bzel.IEN.size());
    vector<array<array<double, dim>, dim>> dN2dx2;
    bzel.Para2Phys(u, v, pt);
    BasisFunction2D(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx, dNdx, detJ);
    for (uint i = 0; i < result_num; i++)
        result[i] = 0.0;

    int count = 0;
    for (uint i = 0; i < bzel.IEN.size(); i++)
    {
        for (uint j = 0; j < state_num; j++)
        {
            result[j] += Nx[i] * val_desire[j][bzel.IEN[i] + time * nPoint];
        }
    }
}