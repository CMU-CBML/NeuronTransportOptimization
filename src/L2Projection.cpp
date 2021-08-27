#include "L2Projection.h"

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
    // cout << "x0: " << x0 << " x1: " << x1 << endl;
    if (x0 >= 0 & x0 <= 0.5 & x1 >= 0 & x1 <= 0.5)
    {
        // result = pow((2.0 * x0 - 1.0), 2.0) * pow((2.0 * x1 - 1.0), 2.0); // test1 and test2
        result = 1.0; // test3 and test4
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

        //DebugMatVector(Stabletmp, e, work_dir + "debug/Stabletmp.txt");

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

void L2Projection::RunL2Projection()
{
    //! L2-projection of desire function
    ReadBezierElementProcess(work_dir);
    nPoint = pts.size();

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
    if (dim == 2)
        BuildLinearSystemProcess2D(pts);
    else
        BuildLinearSystemProcess3D(pts_3d);

    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(Rhs);

    DebugVisualizeVec(Rhs, work_dir + "debug/Rhs_L2.m");
    DebugVisualizeMat(M, work_dir + "debug/M_L2.m");
    //VecAssemblyEnd(Rhs_ini);
    KSPCreate(PETSC_COMM_WORLD, &ksp_L2);
    KSPSetOperators(ksp_L2, M, M);

    KSPGetPC(ksp_L2, &pc_L2);
    PCSetType(pc_L2, PCJACOBI);
    KSPSetTolerances(ksp_L2, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    KSPSolve(ksp_L2, Rhs, X);
    DebugVisualizeVec(X, work_dir + "debug/X_L2.m");
    CollectAndUpdateResult(X);
    //KSPSolve(ksp, Rhs_ini, Y_ini);
}
