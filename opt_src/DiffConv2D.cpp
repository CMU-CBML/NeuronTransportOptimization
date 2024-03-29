#include "DiffConv2D.h"
#include "Utils.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

DiffConv2D::DiffConv2D()
{
    comm = MPI_COMM_WORLD;
    mpiErr = MPI_Comm_rank(comm, &comRank);
    mpiErr = MPI_Comm_size(comm, &comSize);
    nProcess = comSize;
}

DiffConv2D::~DiffConv2D()
{
    // if (Y_d != NULL)
    // {
    //     VecDestroy(&Y_d);
    // }
    // if (Y_bc != NULL)
    // {
    //     VecDestroy(&Y_bc);
    // }
    // if (Y_ini != NULL)
    // {
    //     VecDestroy(&Y_ini);
    // }
    // if (U_ini != NULL)
    // {
    //     VecDestroy(&U_ini);
    // }
    // if (L_ini != NULL)
    // {
    //     VecDestroy(&L_ini);
    // }
    // if (Y_k != NULL)
    // {
    //     VecDestroy(&Y_k);
    // }
    // if (U_k != NULL)
    // {
    //     VecDestroy(&U_k);
    // }
    // if (L_k != NULL)
    // {
    //     VecDestroy(&L_k);
    // }
    // if (dX != NULL)
    // {
    //     VecDestroy(&dX);
    // }
    // if (X != NULL)
    // {
    //     VecDestroy(&X);
    // }
    // if (Res_nl != NULL)
    // {
    //     VecDestroy(&Res_nl);
    // }
    // if (Res_Desire != NULL)
    // {
    //     VecDestroy(&Res_Desire);
    // }
    // if (ResVec != NULL)
    // {
    //     VecDestroy(&ResVec);
    // }
    // VecDestroy(&bsub[0]);
    // VecDestroy(&bsub[1]);
    // VecDestroy(&bsub[2]);

    // MatDestroy(&M);
    // MatDestroy(&K);
    // MatDestroy(&Stable);
    // MatDestroy(&Conv);
    // MatDestroy(&Asubmat[0]);
    // MatDestroy(&Asubmat[1]);
    // MatDestroy(&Asubmat[2]);
    // MatDestroy(&Asubmat[3]);
    // MatDestroy(&Asubmat[4]);
    // MatDestroy(&Asubmat[5]);
    // MatDestroy(&Asubmat[6]);
    // MatDestroy(&Asubmat[7]);
    // MatDestroy(&Asubmat[8]);
    // // MatDestroy(PCsubmat);
    // MatDestroy(&P[0]);
    // MatDestroy(&P[1]);

    // // MatDestroy(&PCMat_tmp);
    // // MatDestroy(&PCMat);
    // MatDestroy(&TanMat_tmp);
    // MatDestroy(&TanMat);

    // KSPDestroy(&ksp);
}

void DiffConv2D::BodyForce(double x, double y, double &Fx, double &Fy)
{
    Fx = 0;
    Fy = 0;
}

void DiffConv2D::ReadBezierElementProcess(string fn)
{
    string stmp;
    int itmp, itmp1;
    int npts, neles, nfunctions;

    string fname_cmat = fn + "cmat.txt";

    // cout << fname_cmat <<endl;

    string fname_debug = work_dir + "debug/ele_process_rank" + to_string(comRank) + ".txt";
    ofstream fout_debug;
    fout_debug.open(fname_debug.c_str());
    for (int i = 0; i < ele_process.size(); i++)
        fout_debug << ele_process[i] << endl;

    string fname_debug1 = work_dir + "debug/read_beziermesh_rank" + to_string(comRank) + ".txt";
    ofstream fout_debug1;
    fout_debug1.open(fname_debug1.c_str());

    ifstream fin_cmat;
    fin_cmat.open(fname_cmat);
    if (fin_cmat.is_open())
    {
        fin_cmat >> neles;
        // cout << neles << endl;
        bzmesh_process.resize(ele_process.size());
        fout_debug1 << "Process: " << comRank << " local bzelement #: " << bzmesh_process.size() << endl;

        int add(0);
        // for (int i = 0; i < neles; i++)
        // {
        //     int itmp, itmp1;
        //     fin_cmat >> itmp >> nfunctions >> itmp1;
        //     // if (comRak == 0)
        //     fout_debug1 << "Process: " << comRank << " itmp: " << itmp << " nfunctions: " << nfunctions << endl;
        //     if (find(ele_process.begin(), ele_process.end(), itmp) != ele_process.end())
        //     {
        //         // if (comRank == 0)
        //         fout_debug1 << "Process: " << comRank << " local element id: " << add << " global element id: " << itmp << endl;
        //         bzmesh_process[add].type = itmp1;
        //         bzmesh_process[add].cmat.resize(nfunctions);
        //         bzmesh_process[add].IEN.resize(nfunctions);
        //         for (int j = 0; j < nfunctions; j++)
        //             fin_cmat >> bzmesh_process[add].IEN[j];
        //         for (int j = 0; j < nfunctions; j++)
        //         {
        //             for (int k = 0; k < 16; k++)
        //             {
        //                 fin_cmat >> bzmesh_process[add].cmat[j][k];
        //             }
        //         }
        //         add++;
        //     }
        //     else
        //     {
        //         // if (comRank == 0)
        //         //     cout << "Process: " << comRank << " global element id: " << itmp << endl;
        //         for (int j = 0; j < nfunctions + 2; j++)
        //             getline(fin_cmat, stmp);
        //     }
        // }
        for (int i = 0; i < neles; i++)
        {
            if (find(ele_process.begin(), ele_process.end(), i) != ele_process.end())
            {
                fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
                fout_debug1 << "Process: " << comRank << " itmp: " << itmp << " nfunctions: " << nfunctions << endl;

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
                if (comRank == 1)
                    fout_debug1 << "Process: " << comRank << " local element id: " << add << endl;
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

void DiffConv2D::GaussInfo(int ng)
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

void DiffConv2D::InitializeProblem(const UserSetting2D *ctx)
{

    MPI_Barrier(comm);

    PetscPrintf(PETSC_COMM_WORLD, "Initializing...\n");

    work_dir = ctx->work_dir;

    /*Initialize parameters*/
    GaussInfo(4);
    n_bzmesh = ctx->n_bzmesh;

    // Scale of the problem
    nPoint = ctx->pts.size();

    // constant parameters
    alpha0 = ctx->var[11];
    alpha1 = ctx->var[12];
    alpha2 = ctx->var[13];
    beta1 = ctx->var[14];
    beta2 = ctx->var[15];
    dt = ctx->var[16];
    nTstep = ctx->var[17];
    par = ctx->var; // Dn0, v_plus, v_minus, k+, k-,k'+,k'-

    // state variables
    for (int i = 0; i < state_num; i++)
    {
        State[i].resize(nPoint * nTstep, 0.0);
    }
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {

            if (ctx->bc_flag[i] == -1)
            {
                State[0][i + j * nPoint] = ctx->val_bc[0][i];
                // State[1][i+j*nPoint] = ctx->val_bc[1][i];
                continue;
            }
        }
    }

    // control variables
    for (int i = 0; i < ctrl_num; i++)
    {
        Ctrl[i].resize(nPoint * nTstep, 0.0);
    }

    // penalty variables
    for (int i = 0; i < state_num; i++)
    {
        Lambda[i].resize(nPoint * nTstep, 0.0);
    }

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {

            if (ctx->bc_flag[i] == -1)
            {
                Lambda[0][i + j * nPoint] = 0.0;
                continue;
            }
        }
    }

    // cout << "nPoint: " << nPoint << endl;
    // cout << "nTstep: " << nTstep << endl;

    PetscPrintf(PETSC_COMM_WORLD, "nPoint: %d\n", nPoint);
    PetscPrintf(PETSC_COMM_WORLD, "nTstep: %d\n", nTstep);
    /*Initialize elements assigned to this rank*/
    AssignProcessor(ctx);

    /*Initialize petsc IS, vector, matrix*/

    ierr = MatCreate(PETSC_COMM_WORLD, &M);
    ierr = MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(M, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(M, 64, NULL, 64, NULL);
    ierr = MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(M);
    ierr = PetscObjectSetName((PetscObject)M, "M");

    ierr = MatCreate(PETSC_COMM_WORLD, &K);
    ierr = MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(K, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(K, 64, NULL, 64, NULL);
    ierr = MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(K);
    ierr = PetscObjectSetName((PetscObject)K, "K");

    ierr = MatCreate(PETSC_COMM_WORLD, &P[0]);
    ierr = MatSetSizes(P[0], PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(P[0], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(P[0], 64, NULL, 64, NULL);
    ierr = MatSetOption(P[0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(P[0]);
    ierr = PetscObjectSetName((PetscObject)P[0], "Px");

    ierr = MatCreate(PETSC_COMM_WORLD, &P[1]);
    ierr = MatSetSizes(P[1], PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(P[1], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(P[1], 64, NULL, 64, NULL);
    ierr = MatSetOption(P[1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(P[1]);
    ierr = PetscObjectSetName((PetscObject)P[1], "Py");

    ierr = MatCreate(PETSC_COMM_WORLD, &Conv);
    ierr = MatSetSizes(Conv, PETSC_DECIDE, PETSC_DECIDE, nPoint, nPoint);
    ierr = MatSetType(Conv, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Conv, 64, NULL, 64, NULL);
    ierr = MatSetOption(Conv, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(Conv);
    ierr = PetscObjectSetName((PetscObject)Conv, "Conv");

    ierr = MatCreate(PETSC_COMM_WORLD, &Stable);
    ierr = MatSetSizes(Stable, PETSC_DECIDE, PETSC_DECIDE,
                       nPoint * nTstep * state_num, nPoint * nTstep * state_num);
    ierr = MatSetType(Stable, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Stable, nPoint, NULL, nPoint, NULL);
    ierr = MatSetOption(Stable, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(Stable);
    ierr = PetscObjectSetName((PetscObject)Stable, "MatStable");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[0]);
    ierr = MatSetSizes(Asubmat[0], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[0], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[0], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[0]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[0], "Asub11");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[1]);
    ierr = MatSetSizes(Asubmat[1], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[1], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[1], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[1]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[1], "Asub12");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[2]);
    ierr = MatSetSizes(Asubmat[2], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[2], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[2], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[2]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[2], "Asub13");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[3]);
    ierr = MatSetSizes(Asubmat[3], PETSC_DECIDE, PETSC_DECIDE,
                       ctrl_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[3], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[3], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[3]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[3], "Asub21");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[4]);
    ierr = MatSetSizes(Asubmat[4], PETSC_DECIDE, PETSC_DECIDE,
                       ctrl_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[4], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[4], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[4]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[4], "Asub22");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[5]);
    ierr = MatSetSizes(Asubmat[5], PETSC_DECIDE, PETSC_DECIDE,
                       ctrl_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[5], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[5], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[5]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[5], "Asub23");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[6]);
    ierr = MatSetSizes(Asubmat[6], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[6], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[6], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[6]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[6], "Asub31");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[7]);
    ierr = MatSetSizes(Asubmat[7], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, ctrl_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[7], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[7], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[7]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[7], "Asub32");

    ierr = MatCreate(PETSC_COMM_WORLD, &Asubmat[8]);
    ierr = MatSetSizes(Asubmat[8], PETSC_DECIDE, PETSC_DECIDE,
                       state_num * nPoint * nTstep, state_num * nPoint * nTstep);
    ierr = MatSetType(Asubmat[8], MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(Asubmat[8], nPoint, NULL, nPoint, NULL);
    ierr = MatSetUp(Asubmat[8]);
    ierr = PetscObjectSetName((PetscObject)Asubmat[8], "Asub33");

    for (int i = 0; i < 9; i++)
    {
        ierr = MatZeroEntries(Asubmat[i]);
    }

    ierr = MatCreate(PETSC_COMM_WORLD, &TanMat);
    ierr = MatSetSizes(TanMat, PETSC_DECIDE, PETSC_DECIDE,
                       (2 * state_num + ctrl_num) * nPoint * nTstep,
                       (2 * state_num + ctrl_num) * nPoint * nTstep);
    ierr = MatSetType(TanMat, MATMPIAIJ);
    ierr = MatMPIAIJSetPreallocation(TanMat, nPoint * 3, NULL, nPoint * 3, NULL);
    ierr = MatSetUp(TanMat);
    ierr = PetscObjectSetName((PetscObject)TanMat, "TanMat");

    PetscPrintf(PETSC_COMM_WORLD, "Setup Initial Vector\n");

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &Y_k);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * ctrl_num * nTstep, &U_k);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &L_k);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &Y_d);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &Y_bc);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * (2 * state_num + ctrl_num) * nTstep, &X);

    ierr = VecSet(Y_k, 1.0);
    ierr = VecSet(U_k, 1.0);
    ierr = VecSet(L_k, 1.0);
    ierr = VecSet(Y_d, 0.0);
    ierr = VecSet(Y_bc, 0.0);
    ierr = VecSet(X, 0.0);

    // vector<int> isg_local[3];
    // PetscInt start, end;

    // ierr = VecGetOwnershipRange(X, &start, &end);
    // PetscPrintf(PETSC_COMM_SELF, "XStart: %d ", start);
    // PetscPrintf(PETSC_COMM_SELF, "Process%d, XStart: %d, XEnd: %d\n", comRank, start, end);
    // for (int i = start; i < end; i++)
    // {
    //     if (i < state_num * nPoint * nTstep)
    //     {
    //         isg_local[0].push_back(i);
    //     }
    //     else if (i < (state_num + ctrl_num) * nPoint * nTstep)
    //     {
    //         isg_local[1].push_back(i);
    //     }
    //     else
    //     {
    //         isg_local[2].push_back(i);
    //     }
    // }
    // PetscInt *isg_val0 = isg_local[0].data();
    // PetscInt *isg_val1 = isg_local[1].data();
    // PetscInt *isg_val2 = isg_local[2].data();

    // ISCreateGeneral(PETSC_COMM_WORLD, isg_local[0].size(), isg_val0, PETSC_COPY_VALUES, &isg[0]);
    // ISCreateGeneral(PETSC_COMM_WORLD, isg_local[1].size(), isg_val1, PETSC_COPY_VALUES, &isg[1]);
    // ISCreateGeneral(PETSC_COMM_WORLD, isg_local[2].size(), isg_val2, PETSC_COPY_VALUES, &isg[2]);

    // DebugVisualizeIS(isg[0], work_dir + "debug/isg0.m");
    // DebugVisualizeIS(isg[1], work_dir + "debug/isg1.m");
    // DebugVisualizeIS(isg[2], work_dir + "debug/isg2.m");

    PetscInt *col_yk, *col_yd;
    PetscReal *vals_yk, *vals_yd;

    // ! Setup Desire State Vec
    PetscMalloc1(nPoint * state_num * nTstep, &col_yd);
    PetscMalloc1(nPoint * state_num * nTstep, &vals_yd);

    int count = 0;
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < nPoint; i++)
        {
            for (int k = 0; k < state_num; k++)
            {
                col_yd[count] = i + k * nPoint + j * state_num * nPoint;
                vals_yd[count] = ctx->val_desire[k][i];
                count++;
            }
        }
    }
    VecSetValues(Y_d, nPoint * state_num * nTstep, col_yd, vals_yd,
                 INSERT_VALUES);
    VecAssemblyBegin(Y_d);
    VecAssemblyEnd(Y_d);

    DebugVisualizeVec(Y_d, work_dir + "debug/Y_d.m");

    PetscFree(col_yd);
    PetscFree(vals_yd);
    PetscPrintf(PETSC_COMM_WORLD, "Set Y_d Done!\n");

    // ! Setup Initial State Vec with BC
    count = 0;
    for (int i = 0; i < ctx->bc_flag.size(); i++)
    {
        if (ctx->bc_flag[i] == -1)
            count += 1;
    }
    n_bcval = count;

    PetscPrintf(PETSC_COMM_WORLD, "BC val: %d\n", count);
    PetscMalloc1(count * nTstep, &col_yk);
    PetscMalloc1(count * nTstep, &vals_yk);
    count = 0;
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {
            if (ctx->bc_flag[i] == -1)
            {
                col_yk[count + 0] = i + 0 * nPoint + j * state_num * nPoint;
                vals_yk[count + 0] = ctx->val_bc[0][i];
                // col_yk[count + 1] = i + 1 * nPoint + j * state_num * nPoint;
                // vals_yk[count + 1] = ctx->val_bc[1][i];
                count += 1;
                continue;
            }
        }
    }
    VecSetValues(Y_k, count, col_yk, vals_yk, INSERT_VALUES);
    VecAssemblyBegin(Y_k);
    VecAssemblyEnd(Y_k);
    VecSetValues(Y_bc, count, col_yk, vals_yk, INSERT_VALUES);
    VecAssemblyBegin(Y_bc);
    VecAssemblyEnd(Y_bc);
    VecSetValues(X, count, col_yk, vals_yk, INSERT_VALUES);
    VecAssemblyBegin(X);
    VecAssemblyEnd(X);
    PetscFree(col_yk);
    PetscFree(vals_yk);
    PetscPrintf(PETSC_COMM_WORLD, "Set X and Y_k Done!\n");

    DebugVisualizeVec(X, work_dir + "debug/X.m");
    DebugVisualizeVec(Y_k, work_dir + "debug/Y_k.m");

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &Res_nl);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * state_num * nTstep, &Res_Desire);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * (2 * state_num + ctrl_num) * nTstep, &dX);
    ierr = VecSet(Res_nl, 0.0);
    ierr = VecSet(Res_Desire, 0.0);

    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,
                        nPoint * nTstep * (state_num * 2 + ctrl_num), &ResVec);

    // ierr = VecCreate(PETSC_COMM_WORLD, &ResVec);
    // ierr = VecSetSizes(ResVec, PETSC_DECIDE, nPoint * nTstep * (state_num * 2 +
    // ctrl_num)); ierr = VecSetFromOptions(ResVec);

    PetscPrintf(PETSC_COMM_WORLD, "Initialization Done!\n");
}

void DiffConv2D::AssignProcessor(const UserSetting2D *ctx)
{
    // Assign the partitioned bezier elements to this processor
    for (int i = 0; i < ctx->ele_process[comRank].size(); i++)
        ele_process.push_back(ctx->ele_process[comRank][i]);
}

/// Need to pay attention that basis function right now is computed for iga
/// using Bezier element definition
void DiffConv2D::BasisFunction(double u, double v, int nen,
                               const vector<array<double, phy_dim>> &pt,
                               const vector<array<double, bzpt_num>> &cmat,
                               vector<double> &Nx,
                               vector<array<double, dim>> &dNdx,
                               vector<array<array<double, dim>, dim>> &dN2dx2,
                               double dudx[dim][dim], double &detJ)
{
    double Nu[4] = {(1. - u) * (1. - u) * (1. - u), 3. * (1. - u) * (1. - u) * u,
                    3. * (1. - u) * u * u, u * u * u};
    double Nv[4] = {(1. - v) * (1. - v) * (1. - v), 3. * (1. - v) * (1. - v) * v,
                    3. * (1. - v) * v * v, v * v * v};

    double dNdu[4] = {-3. * (1. - u) * (1. - u), 3. - 12. * u + 9. * u * u,
                      3. * (2. - 3. * u) * u, 3. * u * u};
    double dNdv[4] = {-3. * (1. - v) * (1. - v), 3. - 12. * v + 9. * v * v,
                      3. * (2. - 3. * v) * v, 3. * v * v};

    double dN2du2[4] = {6. * (1. - u), -12. + 18. * u, 6. - 18. * u, 6. * u};
    double dN2dv2[4] = {6. * (1. - v), -12. + 18. * v, 6. - 18. * v, 6. * v};

    double dNdt[bzpt_num][dim];
    double dN2dt2[bzpt_num][dim][dim];
    double Nx_bz[bzpt_num];
    double dNdx_bz[bzpt_num][dim];
    double dN2dx2_bz[bzpt_num][dim][dim];

    Nx.clear();
    dNdx.clear();
    dN2dx2.clear();
    Nx.resize(nen, 0);
    dNdx.resize(nen, {0});
    dN2dx2.resize(nen, {{0}});

    int i, j, k, a, b, c, loc;
    loc = 0;

    for (j = 0; j < 4; j++)
    {
        for (k = 0; k < 4; k++)
        {
            Nx_bz[loc] = Nu[k] * Nv[j];
            dNdt[loc][0] = dNdu[k] * Nv[j];
            dNdt[loc][1] = Nu[k] * dNdv[j];
            dN2dt2[loc][0][0] = dN2du2[k] * Nv[j];
            dN2dt2[loc][0][1] = dNdu[k] * dNdv[j];
            dN2dt2[loc][1][0] = dNdu[k] * dNdv[j];
            dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j];
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

    // 1st derivatives
    for (i = 0; i < bzpt_num; i++)
    {
        dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
        dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
    }
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            dudx[i][j] = dtdx[i][j];

    detJ = MatrixDet(dxdt);
    detJ = 0.25 * detJ;

    // 2nd derivatives
    double dx2dt2[2][4] = {{0}};
    double dt2dx2[2][4] = {{0}};
    for (int l = 0; l < 2; l++)
        for (loc = 0; loc < bzpt_num; loc++)
            for (a = 0; a < 2; a++)
                for (b = 0; b < 2; b++)
                    dx2dt2[l][2 * a + b] += pt[loc][l] * dN2dt2[loc][a][b];

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                for (a = 0; a < 2; a++)
                    for (b = 0; b < 2; b++)
                        for (c = 0; c < 2; c++)
                            dt2dx2[c][2 * i + j] -=
                                dx2dt2[k][2 * a + b] * dtdx[a][i] * dtdx[b][j] * dtdx[c][k];

    for (loc = 0; loc < bzpt_num; loc++)
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                dN2dx2_bz[loc][i][j] = 0.;

    for (loc = 0; loc < bzpt_num; loc++)
    {
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                for (a = 0; a < 2; a++)
                {
                    for (b = 0; b < 2; b++)
                    {
                        dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i] * dtdx[b][j];
                    }
                    dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][2 * i + j];
                }
            }
        }
    }

    for (i = 0; i < nen; i++)
    {
        for (j = 0; j < bzpt_num; j++)
        {
            Nx[i] += cmat[i][j] * Nx_bz[j];
            for (int m = 0; m < 2; m++)
            {
                dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
                for (int n = 0; n < 2; n++)
                {
                    dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
                }
            }
        }
    }
}

void DiffConv2D::BasisFunctionCoarseSpace(
    double u, double v, int nen, const vector<array<double, phy_dim>> &pt,
    const vector<array<double, bzpt_num>> &cmat, vector<double> &Nx,
    vector<array<double, dim>> &dNdx, double &detJ)
{
    double Nu[2] = {(1. - u), u};
    double Nv[2] = {(1. - v), v};

    double dNdu[2] = {-1., 1.};
    double dNdv[2] = {-1., 1.};

    double dNdt[4][dim];

    Nx.clear();
    dNdx.clear();
    Nx.resize(nen, 0);
    dNdx.resize(nen, {0});

    int i, j, k, a, b, c, loc;
    loc = 0;

    for (k = 0; k < 2; k++)
    {
        Nx[loc] = Nu[k] * Nv[0];
        dNdt[loc][0] = dNdu[k] * Nv[0];
        dNdt[loc][1] = Nu[k] * dNdv[0];
        loc++;
    }
    for (k = 0; k < 2; k++)
    {
        Nx[loc] = Nu[1 - k] * Nv[1];
        dNdt[loc][0] = dNdu[1 - k] * Nv[1];
        dNdt[loc][1] = Nu[1 - k] * dNdv[1];
        loc++;
    }

    double dxdt[2][2] = {{0}};
    int id_linear[4] = {0, 3, 15, 12};
    for (loc = 0; loc < 4; loc++)
        for (a = 0; a < 2; a++)
            for (b = 0; b < 2; b++)
                dxdt[a][b] += pt[id_linear[loc]][a] * dNdt[loc][b];

    // for (loc = 0; loc < 4;loc++)
    // {
    // 	cout << "dNdt: " << loc << " " << dNdt[loc][0] << " " << dNdt[loc][1] <<
    // "\n"; 	cout << "pt: " << loc << " " << pt[id_linear[loc]][0] << " " <<
    // pt[id_linear[loc]][1] << "\n";
    // }

    // for (a = 0; a < 2; a++)
    // 	for (b = 0; b < 2; b++)
    // 		cout << "dxdt[a][b]" << dxdt[a][b] << "\n";

    double dtdx[2][2] = {{0}};
    Matrix2DInverse(dxdt, dtdx);

    // 1st derivatives
    for (i = 0; i < 4; i++)
    {
        dNdx[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
        dNdx[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
    }

    detJ = MatrixDet(dxdt);
    detJ = 0.25 * detJ;
}

void DiffConv2D::PointFormValue(vector<double> &Nx, const vector<double> &U,
                                double Value)
{
    Value = 0.;
    for (int j = 0; j < Nx.size(); j++)
        Value += U[j] * Nx[j];
}

void DiffConv2D::PointFormGrad(vector<array<double, dim>> &dNdx,
                               const vector<double> &U, double Value[dim])
{
    for (int j = 0; j < dim; j++)
        Value[j] = 0.;

    for (int j = 0; j < dim; j++)
        for (int k = 0; k < dNdx.size(); k++)
            Value[j] += U[k] * dNdx[k][j];
}

void DiffConv2D::PointFormHess(vector<array<array<double, dim>, dim>> &d2Ndx2,
                               const vector<double> &U,
                               double Value[dim][dim])
{

    for (int j = 0; j < dim; j++)
    {
        for (int k = 0; k < dim; k++)
        {
            Value[j][k] = 0.;
        }
    }

    for (int j = 0; j < dim; j++)
    {
        for (int k = 0; k < dim; k++)
        {
            for (int l = 0; l < d2Ndx2.size(); l++)
            {
                Value[j][k] += U[l] * d2Ndx2[l][j][k];
            }
        }
    }
}

void DiffConv2D::ResnlAssembly(vector<double> Eqtmp, const vector<int> &IEN,
                               Vec &Gqtmp)
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
                vals[i + nen * k + j * nen * state_num] =
                    Eqtmp[i + nen * k + j * nen * state_num];
            }
        }
    }
    VecSetValues(Gqtmp, nen * state_num * nTstep, rows, vals, ADD_VALUES);

    PetscFree(rows);
    PetscFree(vals);

    delete rows, vals;
}

void DiffConv2D::StableMatAssembly(vector<vector<double>> Emat,
                                   const vector<int> &IEN, Mat &Gmat)
{
    int i, j, A, B, m, n;
    int nen = IEN.size();
    int row_start, row_end, row_now;
    int add = 0;

    PetscInt *nodeList;
    PetscReal *tmpGmat;

    PetscMalloc1(nen, &nodeList);
    PetscMalloc1(nen, &tmpGmat);

    // PetscInt *nodeList = new PetscInt[nen * nTstep];
    // PetscReal *tmpGmat = new PetscReal[nen * nen * nTstep * nTstep ];

    for (i = 0; i < state_num; i++)
    {
        for (j = 0; j < nTstep; j++)
        {
            for (m = 0; m < nen; m++)
            {
                A = IEN[m] + i * nPoint + j * nPoint * state_num;
                for (n = 0; n < nen; n++)
                {
                    B = IEN[n] + i * nPoint + j * nPoint * state_num;
                    nodeList[n] = B;
                    tmpGmat[n] = Emat[m + j * nen][n];
                }
                MatSetValues(Gmat, 1, &A, nen, nodeList, tmpGmat, ADD_VALUES);
            }
        }
    }

    PetscFree(nodeList);
    PetscFree(tmpGmat);

    delete nodeList;
    delete tmpGmat;
}

void DiffConv2D::MatrixAssembly(vector<vector<double>> Emat,
                                const vector<int> &IEN, Mat &Gmat)
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

    MatSetValues(Gmat, IEN.size(), nodeList, IEN.size(), nodeList, tmpGmat,
                 ADD_VALUES);
    delete nodeList;
    delete tmpGmat;
}

void DiffConv2D::ComputeMassMatrix(vector<double> &Nx, const double detJ,
                                   vector<vector<double>> &MassMat)
{
    int a, b, nen = Nx.size();
    for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
            MassMat[a][b] += Nx[a] * Nx[b] * detJ;
}

void DiffConv2D::ComputeStiffMatrix(vector<array<double, dim>> &dNdx,
                                    const double detJ,
                                    vector<vector<double>> &StiffMat)
{
    int a, b, c, nen = dNdx.size();
    double tmp = 0;
    for (a = 0; a < nen; a++)
    {
        for (b = 0; b < nen; b++)
        {
            tmp = 0;
            for (c = 0; c < dim; c++)
            {
                tmp += dNdx[a][c] * dNdx[b][c] * detJ;
            }
            StiffMat[a][b] += tmp;
        }
    }
}

void DiffConv2D::ComputeConvectMatrix(vector<double> &v, vector<double> &Nx,
                                      vector<array<double, dim>> &dNdx,
                                      const double detJ,
                                      vector<vector<double>> &ConvectMat)
{
    int a, b, nen = Nx.size();
    double tmp = 0;
    for (a = 0; a < nen; a++)
    {
        for (b = 0; b < nen; b++)
        {
            tmp = 0;
            for (int c = 0; c < dim; c++)
            {
                tmp += Nx[a] * dNdx[b][c] * v[c] * detJ;
            }
            ConvectMat[a][b] += tmp;
        }
    }
}

void DiffConv2D::ComputeParMatrix(vector<double> &Nx,
                                  vector<array<double, dim>> &dNdx,
                                  const double detJ, int dir,
                                  vector<vector<double>> &ParMat)
{
    int a, b, nen = Nx.size();
    for (a = 0; a < nen; a++)
        for (b = 0; b < nen; b++)
            ParMat[a][b] += Nx[a] * dNdx[b][dir] * detJ;
}

void DiffConv2D::LocalL2Projection(int e, vector<double> w,
                                   vector<vector<double>> &u_proj)
{

    int nen, A;
    nen = bzmesh_process[e].IEN.size();
    int ne_disc = 4;

    u_proj.clear();
    u_proj.resize(nen);
    for (int i = 0; i < nen; i++)
    {
        u_proj[i].resize(ne_disc);
    }

    // * Fine space
    double dudx[dim][dim];
    double detJ;
    vector<double> Nx;
    vector<array<double, 2>> dNdx;
    vector<array<array<double, 2>, 2>> dN2dx2;

    // * Coarse discontinuous space

    double detJ_disc, detJ_all(0.0);
    vector<double> Nx_disc;
    vector<array<double, 2>> dNdx_disc;

    vector<vector<double>> Mtmp;
    vector<double> Rhstmp;
    Mat EM;
    Vec Uproj, ERhs;

    MatCreate(PETSC_COMM_SELF, &EM);
    MatSetType(EM, MATAIJ);
    MatSetSizes(EM, PETSC_DECIDE, PETSC_DECIDE, ne_disc, ne_disc);
    MatSetUp(EM);

    VecCreate(PETSC_COMM_SELF, &ERhs);
    VecSetSizes(ERhs, PETSC_DECIDE, ne_disc);
    VecSetUp(ERhs);
    VecDuplicate(ERhs, &Uproj);

    // * Initialize all variables

    Nx.clear();
    Nx.resize(nen, 0);
    dNdx.clear();
    dNdx.resize(nen, {0});
    dN2dx2.clear();
    dN2dx2.resize(nen, {{0}});

    Nx_disc.clear();
    Nx_disc.resize(ne_disc, 0);
    dNdx_disc.clear();
    dNdx_disc.resize(ne_disc, {0});

    for (int idx = 0; idx < nen; idx++)
    {
        // cout << "idx: " << idx << "\n";
        Mtmp.clear();
        Mtmp.resize(ne_disc);
        Rhstmp.clear();
        Rhstmp.resize(ne_disc, 0.0);

        for (int i = 0; i < ne_disc; i++)
        {
            Mtmp[i].resize(ne_disc, 0.);
        }

        // cout << "idx: " << idx << "\n";
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                // * Compute Matrix
                BasisFunctionCoarseSpace(Gpt[i], Gpt[j], 4, bzmesh_process[e].pts,
                                         bzmesh_process[e].cmat, Nx_disc, dNdx_disc,
                                         detJ_disc);
                detJ_disc = wght[i] * wght[j] * detJ_disc;
                ComputeMassMatrix(Nx_disc, detJ_disc, Mtmp);

                // * Compute Rhs Vec
                double utmp = 0.0;
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;

                utmp = w[0] * dNdx[idx][0] + w[1] * dNdx[idx][1];

                for (int k = 0; k < ne_disc; k++)
                {
                    Rhstmp[k] += Nx_disc[k] * utmp * detJ_disc;
                }
            }
        }

        for (int i = 0; i < ne_disc; i++)
        {
            for (int j = 0; j < ne_disc; j++)
            {
                MatSetValue(EM, i, j, Mtmp[i][j], INSERT_VALUES);
            }
            VecSetValue(ERhs, i, Rhstmp[i], INSERT_VALUES);
        }
        MatAssemblyBegin(EM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(EM, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(ERhs);
        VecAssemblyEnd(ERhs);

        KSP ksp_e;
        KSPCreate(PETSC_COMM_SELF, &ksp_e);
        KSPSetOperators(ksp_e, EM, EM);
        KSPSolve(ksp_e, ERhs, Uproj);

        PetscReal *Uproj_array;
        VecGetArray(Uproj, &Uproj_array);
        // cout << "LocalL2Result:\n";

        // // DebugVisualizeMat(EM, work_dir + "debug/EM.m");
        // // DebugVisualizeVec(ERhs, work_dir + "debug/ERhs.m");

        // // cout << "U:" << u[idx] << "\n";
        for (int i = 0; i < ne_disc; i++)
        {
            u_proj[idx][i] = Uproj_array[i];
        }
        VecRestoreArray(Uproj, &Uproj_array);
        KSPDestroy(&ksp_e);
    }

    MatDestroy(&EM);
    VecDestroy(&ERhs);
    VecDestroy(&Uproj);
}

void DiffConv2D::ComputeStableMatrix(
    const vector<double> val_ini[state_num], double h,
    vector<double> wdNdx_proj, vector<double> &Nx,
    vector<array<double, dim>> &dNdx,
    vector<array<array<double, dim>, dim>> &dN2dx2, const double detJ,
    const vector<int> &IEN, vector<vector<double>> &StableMat)
{
    int a, b, nen = IEN.size();
    double Pe, delta;

    // cout << "IEN.size() = " << nen <<endl;

    vector<double> n0_node, nplus_node, nminus_node;
    // vector<double> vplus_node[dim], vminus_node[dim];
    n0_node.clear();
    n0_node.resize(nen, 0.0);
    nplus_node.clear();
    nplus_node.resize(nen, 0.0);

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < nen; i++)
        {
            if (j != 0)
            {
                int A = IEN[i] + (j - 1) * nPoint;
                n0_node[i] = State[0][A];
                nplus_node[i] = State[1][A];
                // nminus_node[i] = State[2][A];
                // vplus_node[0][i] = State[3][A];
                // vplus_node[1][i] = State[4][A];
                // vminus_node[0][i] = State[5][A];
                // vminus_node[1][i] = State[6][A];
            }
            else
            {
                int A = IEN[i];
                n0_node[i] = val_ini[0][A];
                nplus_node[i] = val_ini[1][A];
                // nminus_node[i] = val_ini[2][A];
                // vplus_node[0][i] = val_ini[3][A];
                // vplus_node[1][i] = val_ini[4][A];
                // vminus_node[0][i] = val_ini[5][A];
                // vminus_node[1][i] = val_ini[6][A];
            }
        }

        double n0Val, n0Grad[dim];
        double npVal, npGrad[dim], vpxVal, vpxGrad[dim], vpyVal, vpyGrad[dim];
        double nmVal, nmGrad[dim], vmxVal, vmxGrad[dim], vmyVal, vmyGrad[dim];

        PointFormValue(Nx, n0_node, n0Val);
        PointFormGrad(dNdx, n0_node, n0Grad);

        PointFormValue(Nx, nplus_node, npVal);
        // PointFormValue(Nx, vplus_node[0], vpxVal);
        // PointFormValue(Nx, vplus_node[1], vpyVal);
        PointFormGrad(dNdx, nplus_node, npGrad);

        Pe = h * sqrt(n0Val * n0Val + npVal * npVal) / 2.0 / par[0];
        // * SUPG stabilization term
        if (Pe < 1.0)
            delta = 0.0;
        else
            delta = h / 2.0 / sqrt(n0Val * n0Val + npVal * npVal) * (1.0 - 1.0 / Pe);

        // * LPS stabilization term
        if (Pe < 1.0)
            delta = 0.0;
        else
            delta = h / sqrt(n0Val * n0Val + npVal * npVal);

        // cout << "Pe: " << Pe << " delta: " << delta << endl;
        // getchar();

        for (int i = 0; i < nen; i++)
        {
            for (int k = 0; k < nen; k++)
            {
                // cout << "i: " << i << " k: " << k <<"\n";
                // a = i + nen * 0 + j * nen * state_num;
                // b = k + nen * 0 + j * nen * state_num;
                a = i + nen * 0 + j * nen;
                b = k;
                // * SUPG stabilization term
                StableMat[a][b] += delta * detJ * dt *
                                   ((n0Val * dNdx[i][0] + npVal * dNdx[i][1]) *
                                        (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) -
                                    par[0] * (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) *
                                        (dN2dx2[i][0][0] + dN2dx2[i][1][1]));
                // * LPS stabilization term

                // double conv_i = n0Val * dNdx[i][0] + npVal * dNdx[i][1];
                // double conv_k = n0Val * dNdx[k][0] + npVal * dNdx[k][1];
                // StableMat[a][b] += delta * detJ * dt * (conv_i - wdNdx_proj[i + j *
                // nen] ) * (conv_k - wdNdx_proj[k + j * nen]);

                // a = i + nen * 1 + j * nen * state_num;
                // b = k + nen * 0 + j * nen * state_num;
                // StableMat[a][b] += delta * detJ * dt * ((n0Val * dNdx[i][0] + npVal *
                // dNdx[i][1]) * (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) - par[0] *
                // (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) * (dN2dx2[i][0][0] +
                // dN2dx2[i][1][1]));
            }
        }
    }

    // for (a = 0; a < nen; a++)
    // 	for (b = 0; b < nen; b++)
    // 		StableMat[a][b] += Nx[a] * dNdx[b][dir] * detJ;
}

void DiffConv2D::ComputeStableMatrix_Convection(
    const vector<double> w, double h, vector<double> wdNdx_proj,
    vector<double> &Nx, vector<array<double, dim>> &dNdx,
    vector<array<array<double, dim>, dim>> &dN2dx2, const double detJ,
    const vector<int> &IEN, vector<vector<double>> &StableMat)
{
    int a, b, nen = IEN.size();
    double Pe, delta;
    double scale;

    // cout << "IEN.size() = " << nen <<endl;

    vector<double> n0_node, nplus_node, nminus_node;
    // vector<double> vplus_node[dim], vminus_node[dim];
    n0_node.clear();
    n0_node.resize(nen, 0.0);
    nplus_node.clear();
    nplus_node.resize(nen, 0.0);

    for (int j = 0; j < nTstep; j++)
    {

        // Pe = h * sqrt(w[0] * w[0] + w[0]* w[1]) / 2.0 / par[0];
        // * SUPG stabilization term
        // if (Pe < 1.0)
        // 	delta = 0.0;
        // else
        // 	delta = h / 2.0 / sqrt(w[0] * w[0] + w[0] * w[1]) * (1.0 - 1.0 / Pe);

        Pe = h * sqrt(w[0] * w[0] + w[1] * w[1]) / par[0];
        // * LPS stabilization term
        if (Pe < 1.0)
            delta = 0.0;
        else
            delta = h / sqrt(w[0] * w[0] + w[1] * w[1]);

        // cout << "Pe: " << Pe << " delta: " << delta << endl;
        // getchar();

        for (int i = 0; i < nen; i++)
        {
            for (int k = 0; k < nen; k++)
            {
                // cout << "i: " << i << " k: " << k <<"\n";
                // a = i + nen * 0 + j * nen * state_num;
                // b = k + nen * 0 + j * nen * state_num;
                a = i + nen * 0 + j * nen;
                b = k;
                // * SUPG stabilization term
                // StableMat[a][b] += delta * detJ * ((n0Val * dNdx[i][0] + npVal *
                // dNdx[i][1]) * (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) - par[0] *
                // (n0Val * dNdx[k][0] + npVal * dNdx[k][1]) * (dN2dx2[i][0][0] +
                // dN2dx2[i][1][1]));
                // * LPS stabilization term

                if (time_int != 0)
                    scale = dt;
                else
                    scale = 1.0;

                double conv_i = w[0] * dNdx[i][0] + w[1] * dNdx[i][1];
                double conv_k = w[0] * dNdx[k][0] + w[1] * dNdx[k][1];
                StableMat[a][b] += delta * detJ * scale *
                                   (conv_i - wdNdx_proj[i + j * nen]) *
                                   (conv_k - wdNdx_proj[k + j * nen]);

                // a = i + nen * 1 + j * nen * state_num;
                // b = k + nen * 0 + j * nen * state_num;
                // StableMat[a][b] += delta * detJ * ((w[0] * dNdx[i][0] + w[1] *
                // dNdx[i][1]) * (w[0] * dNdx[k][0] + w[1] * dNdx[k][1]) - par[0] *
                // (w[0] * dNdx[k][0] + npVal * dNdx[k][1]) * (dN2dx2[i][0][0] +
                // dN2dx2[i][1][1]));
            }
        }
    }

    // for (a = 0; a < nen; a++)
    // 	for (b = 0; b < nen; b++)
    // 		StableMat[a][b] += Nx[a] * dNdx[b][dir] * detJ;
}

void DiffConv2D::ComputeStableMatrix_Convection2(
    const vector<double> w, double h, vector<vector<double>> wdNdx_proj,
    vector<double> &Nx, vector<array<double, dim>> &dNdx,
    vector<double> &Nx_disc, vector<array<double, dim>> &dNdx_disc,
    const double detJ, const double detJ_disc, const vector<int> &IEN,
    vector<vector<double>> &StableMat)
{
    int a, b, nen = IEN.size();
    double Pe, delta;
    double scale;

    for (int j = 0; j < nTstep; j++)
    {

        Pe = h * sqrt(w[0] * w[0] + w[1] * w[1]) / par[0] / 2.0;
        // * LPS stabilization term
        if (Pe < 1.0)
        {
            delta = 0.0;
        }
        else
        {
            // delta = h * h ;
            // delta = h * sqrt(w[0] * w[0] + w[1] * w[1]);
            delta = min(h / sqrt(w[0] * w[0] + w[1] * w[1]), h * h / par[0]);
        }
        // cout << "delta: " << delta << "\n";
        // delta = min(h / sqrt(w[0] * w[0] + w[1] * w[1]), h * h / par[0]);
        delta *= 10;

        for (int i = 0; i < nen; i++)
        {
            for (int k = 0; k < nen; k++)
            {
                a = i + nen * 0 + j * nen;
                b = k;
                // * LPS stabilization term

                if (time_int != 0)
                    scale = dt;
                else
                    scale = 1.0;

                double conv_i = w[0] * dNdx[i][0] + w[1] * dNdx[i][1];
                double conv_k = w[0] * dNdx[k][0] + w[1] * dNdx[k][1];
                double conv_proj_i = 0.0, conv_proj_k = 0.0;
                for (int idx = 0; idx < wdNdx_proj[i].size(); idx++)
                {
                    conv_proj_i += wdNdx_proj[i][idx] * Nx_disc[idx];
                    conv_proj_k += wdNdx_proj[k][idx] * Nx_disc[idx];
                }
                // cout << "conv_i - conv_proj_i " << conv_i - conv_proj_i << "\n";
                StableMat[a][b] += delta * detJ * scale * (conv_i - conv_proj_i) *
                                   (conv_k - conv_proj_k);
            }
        }
    }
}

void DiffConv2D::ComputeResVector(const vector<int> bc_flag,
                                  const vector<double> w,
                                  const vector<double> val_ini[state_num],
                                  vector<double> &Nx,
                                  vector<array<double, dim>> &dNdx,
                                  const vector<int> &IEN, const double detJ,
                                  vector<double> &qtmp)
{
    int a, b, nen = IEN.size();
    // PetscInt *rows;
    // PetscReal *vals;

    // PetscMalloc1(IEN.size()*state_num, &rows);
    // PetscMalloc1(IEN.size()*state_num, &vals);

    vector<double> n0_node, nplus_node, nminus_node;
    // vector<double> vplus_node[dim], vminus_node[dim];
    n0_node.clear();
    n0_node.resize(nen, 0.0);
    // nplus_node.clear(); nplus_node.resize(nen, 0.0);
    // nminus_node.clear(); nminus_node.resize(nen, 0);
    // for (int i = 0; i < dim; i++)
    // {
    // 	vplus_node[i].clear();
    // 	vplus_node[i].resize(nen, 0);
    // 	vminus_node[i].clear();
    // 	vminus_node[i].resize(nen, 0);
    // }

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < nen; i++)
        {
            if (j != 0 || time_int == 0)
            {
                int A = IEN[i] + (j - 1) * nPoint;
                n0_node[i] = State[0][A];
            }
            else
            {
                int A = IEN[i];
                n0_node[i] = val_ini[0][A];
            }
        }

        double n0Val, n0Grad[dim];
        double npVal, npGrad[dim], vpxVal, vpxGrad[dim], vpyVal, vpyGrad[dim];
        double nmVal, nmGrad[dim], vmxVal, vmxGrad[dim], vmyVal, vmyGrad[dim];

        PointFormValue(Nx, n0_node, n0Val);
        PointFormGrad(dNdx, n0_node, n0Grad);

        for (int i = 0; i < nen; i++)
        {
            // qtmp[i + nen * 0] += Nx[i] * n0Val * detJ;

            // ! Convection diffusion
            qtmp[i + nen * 0] += 0.0;
            // qtmp[i + nen * 1] += Nx[i] * npVal * detJ;
        }
    }

    // ! Apply BC while assembly vec
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < nen; i++)
        {
            if (bc_flag[IEN[i]] == -1)
            {
                for (int k = 0; k < nen; k++)
                {
                    if (k != i)
                    {
                        // ! Steady Convection diffusion
                        qtmp[k + nen * 0 + j * nen * state_num] -=
                            State[0][IEN[i]] *
                            (par[0] * (dNdx[k][0] * dNdx[i][0] + dNdx[k][1] * dNdx[i][1]) +
                             (w[0] * dNdx[k][0] + w[1] * dNdx[k][1])) *
                            detJ;
                    }
                }
            }
        }
    }
    // PetscFree(rows);
    // PetscFree(vals);

    // delete rows, vals;
}

void DiffConv2D::GetMatrixPosition(int row, int n_var, int &i_point, int &i_var,
                                   int &i_tstep)
{
    i_tstep = row / (n_var * nPoint);
    i_var = row % (n_var * nPoint) / nPoint;
    i_point = row % (n_var * nPoint) % nPoint;
}

void DiffConv2D::BuildLinearSystemProcessUnitMatrix(
    const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num],
    const vector<double> val_ini[state_num])
{
    // Build linear system in each process
    int e;
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        double h;

        int nen, A;
        double ux_bc, uy_bc, uz_bc, p_bc;

        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;
        vector<vector<double>> Mtmp, Ktmp, Pxtmp, Pytmp, Convtmp, Stabletmp;

        // ! LPS method
        int ne_disc = 4;
        vector<double> wdNdx;
        vector<double> wdNdx_proj;
        vector<vector<double>> wdNdx_proj2;
        double detJ_disc, detJ_all(0.0);
        vector<double> Nx_disc;
        vector<array<double, 2>> dNdx_disc;

        vector<double> n0_node, nplus_node, nminus_node;
        vector<double> vplus_node[dim], vminus_node[dim];
        vector<double> fplus_node[dim], fminus_node[dim];

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0);
        dNdx.clear();
        dNdx.resize(nen, {0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0}});
        wdNdx.clear();
        wdNdx.resize(nen * nTstep, 0.0);

        wdNdx_proj.clear();
        wdNdx_proj.resize(ne_disc * nTstep, 0.0);
        Nx_disc.clear();
        Nx_disc.resize(ne_disc, 0);
        dNdx_disc.clear();
        dNdx_disc.resize(ne_disc, {0});

        // n0_node.clear(); n0_node.resize(nen, 0);
        // nplus_node.clear(); nplus_node.resize(nen, 0);
        // nminus_node.clear(); nminus_node.resize(nen, 0);
        // for(int i=0;i<dim;i++)
        // {
        // 	vplus_node[i].clear(); vplus_node[i].resize(nen, 0);
        // 	vminus_node[i].clear(); vminus_node[i].resize(nen, 0);
        // 	fplus_node[i].clear(); fplus_node[i].resize(nen, 0);
        // 	fminus_node[i].clear(); fminus_node[i].resize(nen, 0);
        // }

        Mtmp.clear();
        Mtmp.resize(nen);
        Ktmp.clear();
        Ktmp.resize(nen);
        Pxtmp.clear();
        Pxtmp.resize(nen);
        Pytmp.clear();
        Pytmp.resize(nen);
        Convtmp.clear();
        Convtmp.resize(nen);
        Stabletmp.clear();
        Stabletmp.resize(nen * nTstep);

        for (int i = 0; i < nen; i++)
        {
            Mtmp[i].resize(nen, 0.);
            Ktmp[i].resize(nen, 0.);
            Pxtmp[i].resize(nen, 0.);
            Pytmp[i].resize(nen, 0.);
            Convtmp[i].resize(nen, 0.);
        }
        for (int i = 0; i < nen * nTstep; i++)
        {
            Stabletmp[i].resize(nen, 0.);
        }
        // getchar();

        // ! Convection velocity
        double theta = PI / 4 * 1.0 / 3.0;
        // theta = 2.4;
        double w_norm = 1.0 * 0.0;
        vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};

        LocalL2Projection(e, w, wdNdx_proj2);
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Start computing!\n";
        int debug_count = 0;
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                ComputeMassMatrix(Nx, detJ, Mtmp);
                ComputeStiffMatrix(dNdx, detJ, Ktmp);
                ComputeConvectMatrix(w, Nx, dNdx, detJ, Convtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 0, Pxtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 1, Pytmp);
                // ComputeStableMatrix(val_ini, h, wdNdx_proj, Nx, dNdx, dN2dx2, detJ,
                // bzmesh_process[e].IEN, Stabletmp);
                // ComputeStableMatrix_Convection(w, h, wdNdx_proj, Nx, dNdx, dN2dx2,
                // detJ, bzmesh_process[e].IEN, Stabletmp);

                // for (int k = 0; k < nen;k++)
                // {
                // 	wdNdx[k] = w[0] * dNdx[k][0] + w[1] * dNdx[k][1];
                // }

                BasisFunctionCoarseSpace(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                                         bzmesh_process[e].cmat, Nx_disc, dNdx_disc,
                                         detJ_disc);
                detJ_disc = wght[i] * wght[j] * detJ_disc;
                h = 1.0 / 10.0;
                ComputeStableMatrix_Convection2(w, h, wdNdx_proj2, Nx, dNdx, Nx_disc,
                                                dNdx_disc, detJ, detJ_disc,
                                                bzmesh_process[e].IEN, Stabletmp);
                // if (comRank == debug_rank)
                // cout << "Process " << comRank << " Element: " << ele_process[e] << "
                // :Start assemblying!\n"; ComputeResVector(val_ini, Nx, dNdx,
                // bzmesh_process[e].IEN, detJ); cout << "Process:" << comRank << " out
                // of " << nProcess << " Breakpoint2 for element " << ele_process[e] <<
                // " Loop count: " << debug_count << ".\n";

                debug_count++;
            }
        }
        // MPI_Barrier(comm);

        // DebugMatVector(Stabletmp, e, work_dir + "debug/Stabletmp.txt");

        // Start element matrix assembly
        MatrixAssembly(Mtmp, bzmesh_process[e].IEN, M);
        MatrixAssembly(Ktmp, bzmesh_process[e].IEN, K);
        MatrixAssembly(Pxtmp, bzmesh_process[e].IEN, P[0]);
        MatrixAssembly(Pytmp, bzmesh_process[e].IEN, P[1]);
        MatrixAssembly(Convtmp, bzmesh_process[e].IEN, Conv);
        StableMatAssembly(Stabletmp, bzmesh_process[e].IEN, Stable);
        // getchar();
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Finish assembly!\n";
    }
    cout << "Process " << comRank << " :complete build matrix and vector!\n";

    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(P[0], MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(P[1], MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Conv, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Stable, MAT_FINAL_ASSEMBLY);
    MPI_Barrier(comm);
    // VecAssemblyBegin(Res_nl);
}

void DiffConv2D::BuildLinearSystemProcessAsubmat(
    const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num],
    const vector<double> val_ini[state_num], const bool flag_glb)
{
    // Build linear system in each process
    int e;
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        double h;

        int nen, A;
        double ux_bc, uy_bc, uz_bc, p_bc;

        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;
        vector<vector<double>> Mtmp, Ktmp, Pxtmp, Pytmp, Convtmp, Stabletmp;

        // ! LPS method
        int ne_disc = 4;
        vector<double> wdNdx;
        vector<double> wdNdx_proj;
        vector<vector<double>> wdNdx_proj2;
        double detJ_disc, detJ_all(0.0);
        vector<double> Nx_disc;
        vector<array<double, 2>> dNdx_disc;

        vector<double> n0_node, nplus_node, nminus_node;
        vector<double> vplus_node[dim], vminus_node[dim];
        vector<double> fplus_node[dim], fminus_node[dim];

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0);
        dNdx.clear();
        dNdx.resize(nen, {0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0}});
        wdNdx.clear();
        wdNdx.resize(nen * nTstep, 0.0);

        wdNdx_proj.clear();
        wdNdx_proj.resize(ne_disc * nTstep, 0.0);
        Nx_disc.clear();
        Nx_disc.resize(ne_disc, 0);
        dNdx_disc.clear();
        dNdx_disc.resize(ne_disc, {0});

        // n0_node.clear(); n0_node.resize(nen, 0);
        // nplus_node.clear(); nplus_node.resize(nen, 0);
        // nminus_node.clear(); nminus_node.resize(nen, 0);
        // for(int i=0;i<dim;i++)
        // {
        // 	vplus_node[i].clear(); vplus_node[i].resize(nen, 0);
        // 	vminus_node[i].clear(); vminus_node[i].resize(nen, 0);
        // 	fplus_node[i].clear(); fplus_node[i].resize(nen, 0);
        // 	fminus_node[i].clear(); fminus_node[i].resize(nen, 0);
        // }

        Mtmp.clear();
        Mtmp.resize(nen);
        Ktmp.clear();
        Ktmp.resize(nen);
        Pxtmp.clear();
        Pxtmp.resize(nen);
        Pytmp.clear();
        Pytmp.resize(nen);
        Convtmp.clear();
        Convtmp.resize(nen);
        Stabletmp.clear();
        Stabletmp.resize(nen * nTstep);

        for (int i = 0; i < nen; i++)
        {
            Mtmp[i].resize(nen, 0.);
            Ktmp[i].resize(nen, 0.);
            Pxtmp[i].resize(nen, 0.);
            Pytmp[i].resize(nen, 0.);
            Convtmp[i].resize(nen, 0.);
        }
        for (int i = 0; i < nen * nTstep; i++)
        {
            Stabletmp[i].resize(nen, 0.);
        }
        // getchar();

        // ! Convection velocity
        double theta = PI / 4 * 1.0 / 3.0;
        // theta = 2.4;
        double w_norm = 1.0 * 0.0;
        vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};

        LocalL2Projection(e, w, wdNdx_proj2);
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Start computing!\n";
        int debug_count = 0;
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                ComputeMassMatrix(Nx, detJ, Mtmp);
                ComputeStiffMatrix(dNdx, detJ, Ktmp);
                ComputeConvectMatrix(w, Nx, dNdx, detJ, Convtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 0, Pxtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 1, Pytmp);
                // ComputeStableMatrix(val_ini, h, wdNdx_proj, Nx, dNdx, dN2dx2, detJ,
                // bzmesh_process[e].IEN, Stabletmp);
                // ComputeStableMatrix_Convection(w, h, wdNdx_proj, Nx, dNdx, dN2dx2,
                // detJ, bzmesh_process[e].IEN, Stabletmp);

                // for (int k = 0; k < nen;k++)
                // {
                // 	wdNdx[k] = w[0] * dNdx[k][0] + w[1] * dNdx[k][1];
                // }

                BasisFunctionCoarseSpace(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                                         bzmesh_process[e].cmat, Nx_disc, dNdx_disc,
                                         detJ_disc);
                detJ_disc = wght[i] * wght[j] * detJ_disc;
                h = 1.0 / 10.0;
                ComputeStableMatrix_Convection2(w, h, wdNdx_proj2, Nx, dNdx, Nx_disc,
                                                dNdx_disc, detJ, detJ_disc,
                                                bzmesh_process[e].IEN, Stabletmp);
                // if (comRank == debug_rank)
                // cout << "Process " << comRank << " Element: " << ele_process[e] << "
                // :Start assemblying!\n"; ComputeResVector(val_ini, Nx, dNdx,
                // bzmesh_process[e].IEN, detJ); cout << "Process:" << comRank << " out
                // of " << nProcess << " Breakpoint2 for element " << ele_process[e] <<
                // " Loop count: " << debug_count << ".\n";

                debug_count++;
            }
        }
        // MPI_Barrier(comm);

        // DebugMatVector(Stabletmp, e, work_dir + "debug/Stabletmp.txt");

        MatrixAssembly(Mtmp, bzmesh_process[e].IEN, M);

        // Start element matrix assembly

        ElementFormMatrixA11(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA22(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA23(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA31(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);

        if (flag_glb)
        {
            ElementFormMatrixA32(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                                 Pxtmp, Pytmp, Stabletmp);
            ElementFormMatrixA13(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                                 Pxtmp, Pytmp, Stabletmp);
            ElementFormMatrixA33(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                                 Pxtmp, Pytmp, Stabletmp);
        }

        // getchar();
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Finish assembly!\n";
    }
    cout << "Process " << comRank << " :complete build matrix and vector!\n";

    if (!flag_glb)
    {
        MatAssemblyBegin(Asubmat[0], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[4], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[5], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[6], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    }
    else
    {
        MatAssemblyBegin(TanMat, MAT_FINAL_ASSEMBLY);
    }

    MPI_Barrier(comm);
}

void DiffConv2D::BuildLinearSystemProcess(
    const UserSetting2D *user, const vector<Vertex2D> &cpts, const vector<double> val_bc[state_num],
    const vector<double> val_ini[state_num], const bool flag_glb)
{
    // Build linear system in each process
    int e;
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        double h;

        int nen, A;
        double ux_bc, uy_bc, uz_bc, p_bc;

        // ! Basis function
        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;

        // ! Unit matrix computation
        vector<vector<double>> Mtmp, Ktmp, Pxtmp, Pytmp, Convtmp, Stabletmp;

        // ! Rhs Vector computation
        vector<double> vec_desire_tmp, b1tmp, b2tmp, b3tmp;
        double tmp_desire[state_num];

        // ! LPS method
        int ne_disc = 4;
        vector<double> wdNdx;
        vector<double> wdNdx_proj;
        vector<vector<double>> wdNdx_proj2;
        double detJ_disc, detJ_all(0.0);
        vector<double> Nx_disc;
        vector<array<double, 2>> dNdx_disc;

        vector<double> n0_node, nplus_node, nminus_node;
        vector<double> vplus_node[dim], vminus_node[dim];
        vector<double> fplus_node[dim], fminus_node[dim];

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0);
        dNdx.clear();
        dNdx.resize(nen, {0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0}});
        wdNdx.clear();
        wdNdx.resize(nen * nTstep, 0.0);

        wdNdx_proj.clear();
        wdNdx_proj.resize(ne_disc * nTstep, 0.0);
        Nx_disc.clear();
        Nx_disc.resize(ne_disc, 0);
        dNdx_disc.clear();
        dNdx_disc.resize(ne_disc, {0});

        // n0_node.clear(); n0_node.resize(nen, 0);
        // nplus_node.clear(); nplus_node.resize(nen, 0);
        // nminus_node.clear(); nminus_node.resize(nen, 0);
        // for(int i=0;i<dim;i++)
        // {
        // 	vplus_node[i].clear(); vplus_node[i].resize(nen, 0);
        // 	vminus_node[i].clear(); vminus_node[i].resize(nen, 0);
        // 	fplus_node[i].clear(); fplus_node[i].resize(nen, 0);
        // 	fminus_node[i].clear(); fminus_node[i].resize(nen, 0);
        // }

        Mtmp.clear();
        Mtmp.resize(nen);
        Ktmp.clear();
        Ktmp.resize(nen);
        Pxtmp.clear();
        Pxtmp.resize(nen);
        Pytmp.clear();
        Pytmp.resize(nen);
        Convtmp.clear();
        Convtmp.resize(nen);
        Stabletmp.clear();
        Stabletmp.resize(nen * nTstep);

        for (int i = 0; i < nen; i++)
        {
            Mtmp[i].resize(nen, 0.);
            Ktmp[i].resize(nen, 0.);
            Pxtmp[i].resize(nen, 0.);
            Pytmp[i].resize(nen, 0.);
            Convtmp[i].resize(nen, 0.);
        }
        for (int i = 0; i < nen * nTstep; i++)
        {
            Stabletmp[i].resize(nen, 0.);
        }

        vec_desire_tmp.resize(nen * state_num * nTstep, 0.0);
        b1tmp.resize(nen * state_num * nTstep, 0.0);
        b2tmp.resize(nen * ctrl_num * nTstep, 0.0);
        b3tmp.resize(nen * state_num * nTstep, 0.0);
        // getchar();

        // ! Convection velocity
        double theta = PI / 4 * 1.0 / 3.0;
        // theta = 2.4;
        double w_norm = 1.0 * 0.0;
        vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};

        LocalL2Projection(e, w, wdNdx_proj2);
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Start computing!\n";
        int debug_count = 0;
        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                ComputeMassMatrix(Nx, detJ, Mtmp);
                ComputeStiffMatrix(dNdx, detJ, Ktmp);
                ComputeConvectMatrix(w, Nx, dNdx, detJ, Convtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 0, Pxtmp);
                ComputeParMatrix(Nx, dNdx, detJ, 1, Pytmp);
                // ComputeStableMatrix(val_ini, h, wdNdx_proj, Nx, dNdx, dN2dx2, detJ,
                // bzmesh_process[e].IEN, Stabletmp);
                // ComputeStableMatrix_Convection(w, h, wdNdx_proj, Nx, dNdx, dN2dx2,
                // detJ, bzmesh_process[e].IEN, Stabletmp);

                // for (int k = 0; k < nen;k++)
                // {
                // 	wdNdx[k] = w[0] * dNdx[k][0] + w[1] * dNdx[k][1];
                // }

                BasisFunctionCoarseSpace(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                                         bzmesh_process[e].cmat, Nx_disc, dNdx_disc,
                                         detJ_disc);
                detJ_disc = wght[i] * wght[j] * detJ_disc;
                h = 1.0 / 10.0;
                ComputeStableMatrix_Convection2(w, h, wdNdx_proj2, Nx, dNdx, Nx_disc,
                                                dNdx_disc, detJ, detJ_disc,
                                                bzmesh_process[e].IEN, Stabletmp);
                // if (comRank == debug_rank)
                // cout << "Process " << comRank << " Element: " << ele_process[e] << "
                // :Start assemblying!\n"; ComputeResVector(val_ini, Nx, dNdx,
                // bzmesh_process[e].IEN, detJ); cout << "Process:" << comRank << " out
                // of " << nProcess << " Breakpoint2 for element " << ele_process[e] <<
                // " Loop count: " << debug_count << ".\n";

                // Compute Rhs Vector
                double X(0.0), Y(0.0), Z(0.0), T(0.0);
                double pt[2];
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                bzmesh_process[e].Para2Phys(Gpt[i], Gpt[j], pt);
                X = pt[0];
                Y = pt[1];
                for (int tt = 0; tt < nTstep; tt++)
                {
                    T = tt * dt;
                    for (int k = 0; k < nen; k++)
                    {
                        user->DesireStateFunction(X, Y, Z, T, tmp_desire);
                        for (int s = 0; s < state_num; s++)
                            vec_desire_tmp[tt * state_num * nen + s * nen + k] +=
                                tmp_desire[s] * Nx[k] * detJ;
                    }
                }

                debug_count++;
            }
        }
        // MPI_Barrier(comm);

        // DebugMatVector(Stabletmp, e, work_dir + "debug/Stabletmp.txt");

        MatrixAssembly(Mtmp, bzmesh_process[e].IEN, M);

        // Start element matrix assembly

        ElementFormMatrixA11(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA22(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA23(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA31(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA32(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA13(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);
        ElementFormMatrixA33(flag_glb, bzmesh_process[e].IEN, Mtmp, Ktmp, Convtmp,
                             Pxtmp, Pytmp, Stabletmp);

        // getchar();
        // if (comRank == debug_rank)
        // cout << "Process " << comRank << " Element: " << ele_process[e] << "
        // :Finish assembly!\n";
    }
    cout << "Process " << comRank << " :complete build matrix and vector!\n";

    if (!flag_glb)
    {
        MatAssemblyBegin(Asubmat[0], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[4], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[5], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Asubmat[6], MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    }
    else
    {
        MatAssemblyBegin(TanMat, MAT_FINAL_ASSEMBLY);
    }

    MPI_Barrier(comm);
}

void DiffConv2D::BuildResVectorProcess(const vector<Vertex2D> &cpts,
                                       const vector<double> val_bc[state_num],
                                       const vector<double> val_ini[state_num],
                                       const vector<int> bc_flag)
{
    // Build linear system in each process
    int e;
    PetscPrintf(PETSC_COMM_WORLD, "Computing Residual Vector.\n");
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        int nen, A;

        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;

        vector<double> qtmp;

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0.0);
        dNdx.clear();
        dNdx.resize(nen, {0.0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0.0}});
        qtmp.clear();
        qtmp.resize(nen * state_num * nTstep, 0.0);

        // ! Convection velocity
        double theta = PI / 4 * 1.0 / 3.0;
        // theta = 2.4;
        double w_norm = 1.0 * 0.0;
        vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};

        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                ComputeResVector(bc_flag, w, val_ini, Nx, dNdx, bzmesh_process[e].IEN,
                                 detJ, qtmp);
            }
        }

        ResnlAssembly(qtmp, bzmesh_process[e].IEN, Res_nl);
        // getchar();
    }
    cout << "Process " << comRank << " :complete compute residual vector!\n";
    VecAssemblyBegin(Res_nl);
}

void DiffConv2D::ElementFormMatrixA11(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = 0;
        colOffset = 0;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            PetscInt *rowList;
            PetscInt *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add = 0;

            PetscMalloc1(IEN.size(), &rowList);
            PetscMalloc1(IEN.size(), &colList);
            PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);

            switch (time_int)
            {
            case 0:
                // !steady state
                scale = 1.0;
                break;
            case 1:
                // !trapezoidal rule
                if (idx_time == 0 || idx_time == (nTstep - 1))
                    scale = 0.5 * dt;
                else
                    scale = dt;
                break;
            case 2:
                // !rectangle rule
                if (idx_time == (nTstep - 1))
                    scale = 0.0;
                else
                    scale = dt;
                break;
            }

            switch (idx_var)
            {
            case 0:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            case 1:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            }
            if (!flag_glb)
            {
                MatSetValues(Asubmat[0], IEN.size(), rowList, IEN.size(), colList,
                             tmpGmat, ADD_VALUES);
            }
            else
            {
                MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList, tmpGmat,
                             ADD_VALUES);
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}

void DiffConv2D::ElementFormMatrixA12(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable) {}

void DiffConv2D::ElementFormMatrixA13(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = 0;
        colOffset = (state_num + ctrl_num) * nTstep * nPoint;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            PetscInt *rowList, *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add;

            if (time_int != 0)
            {
                if (idx_time > 0)
                {
                    switch (idx_var)
                    {
                    case 0:
                        add = 0;

                        PetscMalloc1(IEN.size() * 2, &rowList);
                        PetscMalloc1(IEN.size(), &colList);
                        PetscMalloc1(IEN.size() * IEN.size() * 2, &tmpGmat);

                        //  ! Unfinished
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            colList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + colOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                rowList[2 * n + 0] = B + idx_var * nPoint +
                                                     (idx_time - 1) * nPoint * state_num +
                                                     rowOffset;
                                rowList[2 * n + 1] = B + idx_var * nPoint +
                                                     idx_time * nPoint * state_num + rowOffset;
                                tmpGmat[m + (2 * n + 0) * IEN.size()] = -EM[m][n];
                                tmpGmat[m + (2 * n + 1) * IEN.size()] =
                                    EM[m][n] + dt * par[0] * EK[m][n] +
                                    EStable[m + idx_time * IEN.size()][n];
                                add += 2;
                            }
                        }

                        break;
                    case 1:
                        add = 0;
                        PetscMalloc1(IEN.size() * 2, &rowList);
                        PetscMalloc1(IEN.size(), &colList);
                        PetscMalloc1(IEN.size() * IEN.size() * 2, &tmpGmat);
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            colList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + colOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                rowList[2 * n + 0] = B + idx_var * nPoint +
                                                     (idx_time - 1) * nPoint * state_num +
                                                     rowOffset;
                                rowList[2 * n + 1] = B + idx_var * nPoint +
                                                     idx_time * nPoint * state_num + rowOffset;
                                tmpGmat[m + (2 * n + 0) * IEN.size()] = -EM[m][n];
                                tmpGmat[m + (2 * n + 1) * IEN.size()] =
                                    EM[m][n] + dt * par[0] * EK[m][n] +
                                    EStable[m + idx_time * IEN.size()][n];
                                add += 2;
                            }
                        }
                        break;
                    }
                    if (!flag_glb)
                    {
                        MatSetValues(Asubmat[6], IEN.size() * 2, rowList, IEN.size(),
                                     colList, tmpGmat, ADD_VALUES);
                    }
                    else
                    {
                        MatSetValues(TanMat, IEN.size() * 2, rowList, IEN.size(), colList,
                                     tmpGmat, ADD_VALUES);
                    }
                }
                else
                {
                    add = 0;

                    PetscMalloc1(IEN.size(), &rowList);
                    PetscMalloc1(IEN.size(), &colList);
                    PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);
                    switch (idx_var)
                    {
                    case 0:
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            colList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + colOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                rowList[n] = B + idx_var * nPoint +
                                             idx_time * nPoint * state_num + rowOffset;
                                tmpGmat[m + n * IEN.size()] =
                                    EM[m][n] + dt * par[0] * EK[m][n] +
                                    EStable[m + idx_time * IEN.size()][n];
                                add++;
                            }
                        }
                        break;
                    case 1:
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            colList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + colOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                rowList[n] = B + idx_var * nPoint +
                                             idx_time * nPoint * state_num + rowOffset;
                                tmpGmat[m + n * IEN.size()] =
                                    EM[m][n] + dt * par[0] * EK[m][n] +
                                    EStable[m + idx_time * IEN.size()][n];
                                add++;
                            }
                        }
                        break;
                    }
                    if (!flag_glb)
                    {
                        MatSetValues(Asubmat[2], IEN.size(), rowList, IEN.size(), colList,
                                     tmpGmat, ADD_VALUES);
                    }
                    else
                    {
                        MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList,
                                     tmpGmat, ADD_VALUES);
                    }
                }
            }
            else
            {
                add = 0;
                PetscMalloc1(IEN.size(), &rowList);
                PetscMalloc1(IEN.size(), &colList);
                PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    colList[m] = A + idx_var * nPoint + colOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        rowList[n] = B + idx_var * nPoint + rowOffset;
                        // ! Pure Diffusion
                        // tmpGmat[add] = par[0] * EK[m][n];
                        // ! Convection
                        tmpGmat[add] = par[0] * EK[m][n] + EConv[m][n] +
                                       EStable[m + idx_time * IEN.size()][n];
                        // tmpGmat[add] = par[0] * EK[m][n] + EConv[m][n];

                        add++;
                    }
                }
                if (!flag_glb)
                {
                    MatSetValues(Asubmat[2], IEN.size(), rowList, IEN.size(), colList,
                                 tmpGmat, ADD_VALUES);
                }
                else
                {
                    MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList,
                                 tmpGmat, ADD_VALUES);
                }
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}
void DiffConv2D::ElementFormMatrixA21(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable) {}
void DiffConv2D::ElementFormMatrixA22(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = state_num * nTstep * nPoint;
        colOffset = state_num * nTstep * nPoint;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < ctrl_num; idx_var++)
        {
            PetscInt *rowList;
            PetscInt *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add;

            PetscMalloc1(IEN.size(), &rowList);
            PetscMalloc1(IEN.size(), &colList);
            PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);

            switch (time_int)
            {
            case 0:
                // !steady state
                scale = 1.0;
                break;
            case 1:
                // !trapezoidal rule
                if (idx_time == 0 || idx_time == (nTstep - 1))
                    scale = 0.5 * dt;
                else
                    scale = dt;
                break;
            case 2:
                // !rectangle rule
                if (idx_time == (nTstep - 1))
                    scale = 0.0;
                else
                    scale = dt;
                break;
            }

            switch (idx_var)
            {
            case 0:
                scale = scale * beta1;
                add = 0;
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] =
                            B + idx_var * nPoint + idx_time * nPoint * ctrl_num + colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            case 1:
                scale = scale * beta1;
                add = 0;
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] =
                            B + idx_var * nPoint + idx_time * nPoint * ctrl_num + colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            }
            if (!flag_glb)
            {
                MatSetValues(Asubmat[4], IEN.size(), rowList, IEN.size(), colList,
                             tmpGmat, ADD_VALUES);
            }
            else
            {
                MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList, tmpGmat,
                             ADD_VALUES);
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}
void DiffConv2D::ElementFormMatrixA23(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = state_num * nTstep * nPoint;
        colOffset = (state_num + ctrl_num) * nTstep * nPoint;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < ctrl_num; idx_var++)
        {
            PetscInt *rowList, *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add = 0;

            PetscMalloc1(IEN.size(), &rowList);
            PetscMalloc1(IEN.size(), &colList);
            PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);

            switch (time_int)
            {
            case 0:
                // !steady state
                scale = -1.0;
                break;
            case 1:
                // !trapezoidal rule
                scale = -dt;
                break;
            case 2:
                // !rectangle rule
                scale = -dt;
                break;
            }

            switch (idx_var)
            {
            case 0:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            case 1:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     colOffset;
                        tmpGmat[add] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            }
            if (!flag_glb)
            {
                MatSetValues(Asubmat[5], IEN.size(), rowList, IEN.size(), colList,
                             tmpGmat, ADD_VALUES);
            }
            else
            {
                MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList, tmpGmat,
                             ADD_VALUES);
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}

void DiffConv2D::ElementFormMatrixA31(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = (state_num + ctrl_num) * nTstep * nPoint;
        colOffset = 0;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            PetscInt *rowList, *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add;

            if (time_int != 0)
            {
                if (idx_time > 0)
                {
                    switch (idx_var)
                    {
                    case 0:
                        add = 0;

                        PetscMalloc1(IEN.size(), &rowList);
                        PetscMalloc1(IEN.size() * 2, &colList);
                        PetscMalloc1(IEN.size() * IEN.size() * 2, &tmpGmat);

                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            rowList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + rowOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                colList[2 * n + 0] = B + idx_var * nPoint +
                                                     (idx_time - 1) * nPoint * state_num +
                                                     colOffset;
                                colList[2 * n + 1] = B + idx_var * nPoint +
                                                     idx_time * nPoint * state_num + colOffset;
                                tmpGmat[add] = -EM[m][n];
                                tmpGmat[add + 1] = EM[m][n] + dt * par[0] * EK[m][n] +
                                                   EStable[m + idx_time * IEN.size()][n];
                                add += 2;
                            }
                        }
                        break;
                    case 1:
                        add = 0;
                        PetscMalloc1(IEN.size(), &rowList);
                        PetscMalloc1(IEN.size() * 2, &colList);
                        PetscMalloc1(IEN.size() * IEN.size() * 2, &tmpGmat);
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            rowList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + rowOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                colList[2 * n + 0] = B + idx_var * nPoint +
                                                     (idx_time - 1) * nPoint * state_num +
                                                     colOffset;
                                colList[2 * n + 1] = B + idx_var * nPoint +
                                                     idx_time * nPoint * state_num + colOffset;
                                tmpGmat[add] = -EM[m][n];
                                tmpGmat[add + 1] = EM[m][n] + dt * par[0] * EK[m][n] +
                                                   EStable[m + idx_time * IEN.size()][n];
                                add += 2;
                            }
                        }
                        break;
                    }
                    if (!flag_glb)
                    {
                        MatSetValues(Asubmat[6], IEN.size(), rowList, IEN.size() * 2,
                                     colList, tmpGmat, ADD_VALUES);
                    }
                    else
                    {
                        MatSetValues(TanMat, IEN.size(), rowList, IEN.size() * 2, colList,
                                     tmpGmat, ADD_VALUES);
                    }
                }
                else
                {
                    add = 0;

                    PetscMalloc1(IEN.size(), &rowList);
                    PetscMalloc1(IEN.size(), &colList);
                    PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);
                    switch (idx_var)
                    {
                    case 0:
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            rowList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + rowOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                colList[n] = B + idx_var * nPoint +
                                             idx_time * nPoint * state_num + colOffset;
                                tmpGmat[add] = EM[m][n] + dt * par[0] * EK[m][n] +
                                               EStable[m + idx_time * IEN.size()][n];
                                add++;
                            }
                        }
                        break;
                    case 1:
                        for (m = 0; m < IEN.size(); m++)
                        {
                            A = IEN[m];
                            rowList[m] = A + idx_var * nPoint +
                                         idx_time * nPoint * state_num + rowOffset;
                            for (n = 0; n < IEN.size(); n++)
                            {
                                B = IEN[n];
                                colList[n] = B + idx_var * nPoint +
                                             idx_time * nPoint * state_num + colOffset;
                                tmpGmat[add] = EM[m][n] + dt * par[0] * EK[m][n] +
                                               EStable[m + idx_time * IEN.size()][n];
                                add++;
                            }
                        }
                        break;
                    }
                    if (!flag_glb)
                    {
                        MatSetValues(Asubmat[6], IEN.size(), rowList, IEN.size(), colList,
                                     tmpGmat, ADD_VALUES);
                    }
                    else
                    {
                        MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList,
                                     tmpGmat, ADD_VALUES);
                    }
                }
            }
            else
            {
                add = 0;
                PetscMalloc1(IEN.size(), &rowList);
                PetscMalloc1(IEN.size(), &colList);
                PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    rowList[m] = A + idx_var * nPoint + rowOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        colList[n] = B + idx_var * nPoint + colOffset;
                        // ! Pure Diffusion
                        // tmpGmat[add] = par[0] * EK[m][n];
                        // ! Convection
                        tmpGmat[add] = par[0] * EK[m][n] + EConv[m][n] +
                                       EStable[m + idx_time * IEN.size()][n];
                        // tmpGmat[add] = par[0] * EK[m][n] + EConv[m][n];

                        add++;
                    }
                }
                if (!flag_glb)
                {
                    MatSetValues(Asubmat[6], IEN.size(), rowList, IEN.size(), colList,
                                 tmpGmat, ADD_VALUES);
                }
                else
                {
                    MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList,
                                 tmpGmat, ADD_VALUES);
                }
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}

void DiffConv2D::ElementFormMatrixA32(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    if (flag_glb)
    {
        rowOffset = (state_num + ctrl_num) * nTstep * nPoint;
        colOffset = state_num * nTstep * nPoint;
    }
    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < ctrl_num; idx_var++)
        {
            PetscInt *rowList, *colList;
            PetscReal *tmpGmat;
            PetscReal scale;
            int add = 0;

            PetscMalloc1(IEN.size(), &rowList);
            PetscMalloc1(IEN.size(), &colList);
            PetscMalloc1(IEN.size() * IEN.size(), &tmpGmat);

            switch (time_int)
            {
            case 0:
                // !steady state
                scale = -1.0;
                break;
            case 1:
                // !trapezoidal rule
                scale = -dt;
                break;
            case 2:
                // !rectangle rule
                scale = -dt;
                break;
            }

            switch (idx_var)
            {
            case 0:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    colList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + colOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        rowList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     rowOffset;
                        tmpGmat[m + n * IEN.size()] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            case 1:
                for (m = 0; m < IEN.size(); m++)
                {
                    A = IEN[m];
                    colList[m] =
                        A + idx_var * nPoint + idx_time * nPoint * ctrl_num + colOffset;
                    for (n = 0; n < IEN.size(); n++)
                    {
                        B = IEN[n];
                        rowList[n] = B + idx_var * nPoint + idx_time * nPoint * state_num +
                                     rowOffset;
                        tmpGmat[m + n * IEN.size()] = scale * EM[m][n];
                        add++;
                    }
                }
                break;
            }
            if (!flag_glb)
            {
                MatSetValues(Asubmat[7], IEN.size(), rowList, IEN.size(), colList,
                             tmpGmat, ADD_VALUES);
            }
            else
            {
                MatSetValues(TanMat, IEN.size(), rowList, IEN.size(), colList, tmpGmat,
                             ADD_VALUES);
            }
            PetscFree(rowList);
            PetscFree(colList);
            PetscFree(tmpGmat);
        }
    }
}

void DiffConv2D::ElementFormMatrixA33(bool flag_glb, const vector<int> &IEN,
                                      vector<vector<double>> EM,
                                      vector<vector<double>> EK,
                                      vector<vector<double>> EConv,
                                      vector<vector<double>> EPx,
                                      vector<vector<double>> EPy,
                                      vector<vector<double>> EStable)
{
    int i, j, A, B, m, n;
    int rowOffset(0), colOffset(0);
    int row, col;
    if (flag_glb)
    {
        rowOffset = (state_num + ctrl_num) * nTstep * nPoint;
        colOffset = (state_num + ctrl_num) * nTstep * nPoint;
    }

    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        // for (int idx_var = 0; idx_var < state_num; idx_var++)
        // {
        for (m = 0; m < IEN.size(); m++)
        {
            A = IEN[m];
            // row = A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
            row = A + idx_time * nPoint * state_num + rowOffset;
            col = row + colOffset;

            if (!flag_glb)
            {
                MatSetValue(Asubmat[8], row, row, 0.0, ADD_VALUES);
            }
            else
            {
                MatSetValue(TanMat, row, row, 0.0, ADD_VALUES);
            }
            // }
        }
    }
}

void DiffConv2D::ElementFormRhsVecb1(bool flag_glb, const vector<int> &IEN,
                                     vector<vector<double>> EM,
                                     vector<vector<double>> EK,
                                     vector<vector<double>> EConv,
                                     vector<vector<double>> EPx,
                                     vector<vector<double>> EPy,
                                     vector<vector<double>> EStable,
                                     vector<double>
                                         Evec_desire)
{
    int i, j, A, B, m, n;
    int nen = IEN.size();
    int rowOffset(0);
    int row, col;
    if (flag_glb)
    {
        rowOffset = 0;
    }

    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            for (m = 0; m < nen; m++)
            {
                A = IEN[m];
                // row = A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
                row = A + idx_time * nPoint * state_num + rowOffset;
                VecSetValue(ResVec, row, Evec_desire[idx_time * state_num * +idx_var * nen + m], ADD_ALL_VALUES);
            }
        }
    }
}

void DiffConv2D::ElementFormRhsVecb2(bool flag_glb, const vector<int> &IEN,
                                     vector<vector<double>> EM,
                                     vector<vector<double>> EK,
                                     vector<vector<double>> EConv,
                                     vector<vector<double>> EPx,
                                     vector<vector<double>> EPy,
                                     vector<vector<double>> EStable,
                                     vector<double>
                                         Evec_desire)
{
    int i, j, A, B, m, n;
    int nen = IEN.size();
    int rowOffset(0);
    int row, col;
    if (flag_glb)
    {
        rowOffset = state_num * nPoint * nTstep;
    }

    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            for (m = 0; m < nen; m++)
            {
                A = IEN[m];
                // row = A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
                row = A + idx_time * nPoint * state_num + rowOffset;
                VecSetValue(ResVec, row, 0.0, ADD_ALL_VALUES);
            }
        }
    }
}

void DiffConv2D::ElementFormRhsVecb3(bool flag_glb, const vector<int> &IEN,
                                     vector<vector<double>> EM,
                                     vector<vector<double>> EK,
                                     vector<vector<double>> EConv,
                                     vector<vector<double>> EPx,
                                     vector<vector<double>> EPy,
                                     vector<vector<double>> EStable,
                                     vector<double>
                                         Evec_desire)
{
    int i, j, A, B, m, n;
    int nen = IEN.size();
    int rowOffset(0);
    int row, col;
    if (flag_glb)
    {
        rowOffset = (state_num + ctrl_num) * nPoint * nTstep;
    }

    for (int idx_time = 0; idx_time < nTstep; idx_time++)
    {
        for (int idx_var = 0; idx_var < state_num; idx_var++)
        {
            for (m = 0; m < nen; m++)
            {
                A = IEN[m];
                // row = A + idx_var * nPoint + idx_time * nPoint * state_num + rowOffset;
                row = A + idx_time * nPoint * state_num + rowOffset;
                for (n = 0; n < nen; n++)
                {
                }
                VecSetValue(ResVec, row, 0.0, ADD_ALL_VALUES);
            }
        }
    }
}

void DiffConv2D::FormMatrixA11(Mat M, Mat K, Mat &A)
{

    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);

    for (row = start; row < end; row++)
    {
        PetscInt ind_point, ind_var, ind_time;
        PetscReal scale;
        PetscInt *col;
        PetscReal *vals;
        PetscInt ncols_m, ncols_k;
        const PetscReal *vals_m, *vals_k;
        const PetscInt *cols_m, *cols_k;
        PetscInt ncols_a;

        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);

        IS is_row, is_col;
        IS *is_row_array, *is_col_array;
        Mat *Msub, *Ksub;

        PetscMalloc1(1, &is_row_array);
        PetscMalloc1(1, &is_col_array);
        PetscMalloc1(1, &Msub);
        PetscMalloc1(1, &Ksub);
        ISCreateStride(PETSC_COMM_SELF, 1, ind_point, 1, &is_row);
        ISCreateStride(PETSC_COMM_SELF, nPoint, 0, 1, &is_col);
        is_row_array[0] = is_row;
        is_col_array[0] = is_col;

        MatCreateSubMatrices(M, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Msub);
        MatCreateSubMatrices(K, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Ksub);
        // // if (row == start && comRank == debug_rank)
        // 	if (comRank == debug_rank)

        // 	{
        // 		DebugVisualizeISSelf(is_row, work_dir + "debug/is_row.m");
        // 		DebugVisualizeISSelf(is_col, work_dir + "debug/is_col.m");
        // 		DebugVisualizeMatSelf(Msub[0], work_dir + "debug/Msub1.m");
        // 		cout << "Process:" << comRank << " out of " << nProcess << "
        // Output debug info.\n";
        // }
        // cout << "Process:" << comRank << " out of " << nProcess << "
        // Breakpoint1.\n";

        MatGetRow(Ksub[0], 0, &ncols_k, &cols_k, &vals_k);
        MatGetRow(Msub[0], 0, &ncols_m, &cols_m, &vals_m);
        // cout << "Process:" << comRank << " out of " << nProcess << "
        // Breakpoint2.\n";

        // MatGetRow(K,ind_point,&ncols_k, &cols_k,&vals_k);
        // MatGetRow(M,ind_point,&ncols_m, &cols_m,&vals_m);
        // if (row == start && comRank == debug_rank)
        // 	if ( comRank == debug_rank)
        // 	{
        // 		cout << "Ncols_k: " << ncols_k << endl;
        // 		for (int i = 0; i < ncols_k; i++)
        // 			cout << "Col: " << cols_k[i] << " Value: " << vals_k[i]
        // << endl; 		cout << "Ncols_m: " << ncols_m << endl;
        // for (int i = 0; i <
        // ncols_m; i++) 			cout << "Col: " << cols_m[i] << " Value: " << vals_m[i]
        // << endl;

        // 		cout << "Output M_row, K_row Done!\n";
        // }
        // getchar();
        switch (ind_var)
        {
        case 0:
            // MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
            ncols_a = ncols_m;
            PetscMalloc1(ncols_a, &col);
            PetscMalloc1(ncols_a, &vals);
            for (int i = 0; i < ncols_a; i++)
            {
                vals[i] = vals_m[i];
                col[i] = cols_m[i] + ind_var * nPoint + ind_time * nPoint * state_num;
            }
            // scale = beta2;
            // scale = 1.0;
            break;
        case 1:
            // MatGetRow(M,ind_point,NULL,NULL,&vals_Extract);
            ncols_a = ncols_m;
            PetscMalloc1(ncols_a, &col);
            PetscMalloc1(ncols_a, &vals);
            for (int i = 0; i < ncols_a; i++)
            {
                vals[i] = vals_m[i];
                col[i] = cols_m[i] + ind_var * nPoint + ind_time * nPoint * state_num;
            }
            // scale = beta2;
            // scale = 1.0;
            break;
        }

        // if(ind_time == 0 || ind_time == (nTstep - 1))
        // 	scale = scale * 0.5 * dt;
        // else
        // 	scale = scale * dt;
        switch (time_int)
        {
        case 0:
            // !steady state
            scale = 1.0;
            break;
        case 1:
            // !trapezoidal rule
            if (ind_time == 0 || ind_time == (nTstep - 1))
                scale = 0.5 * dt;
            else
                scale = dt;
            break;
        case 2:
            // !rectangle rule
            if (ind_time == (nTstep - 1))
                scale = 0.0;
            else
                scale = dt;
            break;
        }

        for (int i = 0; i < ncols_a; i++)
        {
            vals[i] = vals[i] * scale;
        }
        // for(int i = 0; i < nPoint; i++)
        // {
        // 	vals[i] = vals[i]*scale;
        // 	col[i] = i + ind_var * nPoint + ind_time * nPoint * state_num;
        // }

        MatSetValues(A, 1, &row, ncols_a, col, vals, INSERT_VALUES);

        // MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);

        // if(ind_var>=3)
        // MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
        // else if(ind_var>0)
        // MatRestoreRow(K,ind_point,NULL,NULL,&vals_k);

        MatRestoreRow(Msub[0], 0, NULL, NULL, &vals_m);
        // else if(ind_var>0)
        MatRestoreRow(Ksub[0], 0, NULL, NULL, &vals_k);

        PetscFree(is_row_array);
        PetscFree(is_col_array);
        PetscFree(Msub);
        PetscFree(Ksub);
        PetscFree(col);
        PetscFree(vals);
    }

    // ISDestroy(&is_row);
    // ISDestroy(&is_col);

    // PetscFree(is_row_array);
    // PetscFree(is_col_array);
    // PetscFree(Msub);
    // PetscFree(Ksub);
}

void DiffConv2D::FormMatrixA12(Mat M, Mat K, Mat &A)
{
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

// ! Not used. Compute A31 and transpose to get A13
void DiffConv2D::FormMatrixA13(Mat M, Mat K, Mat P[dim], Mat &A)
{
    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);

    for (row = start; row < end; row++)
    {
        PetscInt ind_point, ind_var, ind_time;
        PetscInt ncols_m, ncols_k, ncols_px, ncols_py;
        const PetscReal *vals_m, *vals_k, *vals_px, *vals_py;
        const PetscInt *cols_m, *cols_k, *cols_px, *cols_py;

        PetscInt *col;
        PetscReal *vals;
        PetscReal scale;

        // cout << "Ncols_k: " << ncols_k <<endl;
        // for(int i =0;i<ncols_k;i++)
        // 	cout << "Col: "<< cols_k[i] << " Value: " << vals_k[i] <<endl;
        // cout << "Ncols_m: " << ncols_m <<endl;
        // for(int i =0;i<ncols_m;i++)
        // 	cout << "Col: "<< cols_m[i] << " Value: " << vals_m[i] <<endl;

        // cout << "Output M_row, K_row Done!\n";

        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
        MatGetRow(M, ind_point, &ncols_m, &cols_m, &vals_m);
        MatGetRow(K, ind_point, &ncols_k, &cols_k, &vals_k);
        MatGetRow(P[0], ind_point, &ncols_px, &cols_px, &vals_px);
        MatGetRow(P[1], ind_point, &ncols_py, &cols_py, &vals_py);
        if (ind_time > 0)
        {
            switch (ind_var)
            {
            case 0:
                PetscMalloc1(4 * ncols_m, &col);
                PetscMalloc1(4 * ncols_m, &vals);
                for (int i = 0; i < ncols_m; i++)
                {
                    vals[i + 0 * ncols_m] = -vals_m[i];
                    vals[i + 1 * ncols_m] = vals_m[i] * (1 + dt * par[3] + dt * par[4]);
                    vals[i + 2 * ncols_m] = vals_m[i] * (-dt * par[3]);
                    vals[i + 3 * ncols_m] = vals_m[i] * (-dt * par[4]);
                    col[i + 0 * ncols_m] =
                        cols_m[i] + 0 * nPoint + (ind_time - 1) * nPoint * state_num;
                    col[i + 1 * ncols_m] =
                        cols_m[i] + 0 * nPoint + ind_time * nPoint * state_num;
                    col[i + 2 * ncols_m] =
                        cols_m[i] + 1 * nPoint + ind_time * nPoint * state_num;
                    col[i + 3 * ncols_m] =
                        cols_m[i] + 2 * nPoint + ind_time * nPoint * state_num;
                }
                MatSetValues(A, 1, &row, ncols_m * 4, col, vals, INSERT_VALUES);
                break;
            case 1:
                PetscMalloc1(3 * ncols_m + ncols_px + ncols_py, &col);
                PetscMalloc1(3 * ncols_m + ncols_px + ncols_py, &vals);
                for (int i = 0; i < ncols_m; i++)
                {
                    vals[i + 0 * ncols_m] = -vals_m[i];
                    vals[i + 1 * ncols_m] = vals_m[i] * (-dt * par[5]);
                    vals[i + 2 * ncols_m] = vals_m[i] * (1 + dt * par[5]);
                    col[i + 0 * ncols_m] =
                        cols_m[i] + 1 * nPoint + (ind_time - 1) * nPoint * state_num;
                    col[i + 1 * ncols_m] =
                        cols_m[i] + 0 * nPoint + ind_time * nPoint * state_num;
                    col[i + 2 * ncols_m] =
                        cols_m[i] + 1 * nPoint + ind_time * nPoint * state_num;
                }
                for (int i = 0; i < ncols_px; i++)
                {
                    vals[i + 3 * ncols_m] = vals_px[i] * dt * par[7] * par[7];
                    col[i + 3 * ncols_m] =
                        cols_px[i] + 3 * nPoint + ind_time * nPoint * state_num;
                }
                for (int i = 0; i < ncols_py; i++)
                {
                    vals[i + 3 * ncols_m + ncols_py] = vals_py[i] * dt * par[7] * par[7];
                    col[i + 3 * ncols_m + ncols_py] =
                        cols_py[i] + 4 * nPoint + ind_time * nPoint * state_num;
                }

                // for(int i = 0; i < nPoint; i++)
                // {
                // 	vals[i + 0 * nPoint] = -vals_m[i];
                // 	vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[5]);
                // 	vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[5]);
                // 	vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
                // 	vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
                // 	col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) *
                // nPoint * state_num; 	col[i + 1 * nPoint] = i + ind_var * nPoint +
                // ind_time * nPoint * state_num; 	col[i + 2 * nPoint] = i +
                // (ind_var +
                // 1) * nPoint + ind_time * nPoint * state_num; 	col[i + 3 * nPoint]
                // = i
                // + (ind_var + 3) * nPoint + ind_time * nPoint * state_num; 	col[i +
                // 4 * nPoint] = i + (ind_var + 4) * nPoint + ind_time * nPoint *
                // state_num;
                // }
                MatSetValues(A, 1, &row, (3 * ncols_m + ncols_px + ncols_py), col, vals,
                             INSERT_VALUES);
                break;
            case 2:
                PetscMalloc1(3 * ncols_m + ncols_px + ncols_py, &col);
                PetscMalloc1(3 * ncols_m + ncols_px + ncols_py, &vals);
                for (int i = 0; i < ncols_m; i++)
                {
                    vals[i + 0 * ncols_m] = -vals_m[i];
                    vals[i + 1 * ncols_m] = vals_m[i] * (-dt * par[6]);
                    vals[i + 2 * ncols_m] = vals_m[i] * (1 + dt * par[6]);
                    col[i + 0 * ncols_m] =
                        cols_m[i] + 2 * nPoint + (ind_time - 1) * nPoint * state_num;
                    col[i + 1 * ncols_m] =
                        cols_m[i] + 0 * nPoint + ind_time * nPoint * state_num;
                    col[i + 2 * ncols_m] =
                        cols_m[i] + 2 * nPoint + ind_time * nPoint * state_num;
                }
                for (int i = 0; i < ncols_px; i++)
                {
                    vals[i + 3 * ncols_m] = vals_px[i] * dt * par[7] * par[7];
                    col[i + 3 * ncols_m] =
                        cols_px[i] + 5 * nPoint + ind_time * nPoint * state_num;
                }
                for (int i = 0; i < ncols_py; i++)
                {
                    vals[i + 3 * ncols_m + ncols_py] = vals_py[i] * dt * par[7] * par[7];
                    col[i + 3 * ncols_m + ncols_py] =
                        cols_py[i] + 6 * nPoint + ind_time * nPoint * state_num;
                }
                // PetscMalloc1(5 * nPoint, &col);
                // PetscMalloc1(5 * nPoint, &vals);
                // for(int i = 0; i < nPoint; i++)
                // {
                // 	vals[i + 0 * nPoint] = -vals_m[i];
                // 	vals[i + 1 * nPoint] = vals_m[i] * (-dt * par[6]);
                // 	vals[i + 2 * nPoint] = vals_m[i] * (1 + dt * par[6]);
                // 	vals[i + 3 * nPoint] = vals_px[i] * dt * par[7] * par[7];
                // 	vals[i + 4 * nPoint] = vals_py[i] * dt * par[7] * par[7];
                // 	col[i + 0 * nPoint] = i + ind_var * nPoint + (ind_time - 1) *
                // nPoint * state_num; 	col[i + 1 * nPoint] = i + ind_var * nPoint +
                // ind_time * nPoint * state_num; 	col[i + 2 * nPoint] = i +
                // (ind_var +
                // 1) * nPoint + ind_time * nPoint * state_num; 	col[i + 3 * nPoint]
                // = i
                // + (ind_var + 5) * nPoint + ind_time * nPoint * state_num; 	col[i +
                // 4 * nPoint] = i + (ind_var + 6) * nPoint + ind_time * nPoint *
                // state_num;
                // }
                MatSetValues(A, 1, &row, (3 * ncols_m + ncols_px + ncols_py), col, vals,
                             INSERT_VALUES);
                // PetscFree(col);
                // PetscFree(vals);
                break;
            default:

                PetscMalloc1(2 * ncols_m, &col);
                PetscMalloc1(2 * ncols_m, &vals);
                for (int i = 0; i < ncols_m; i++)
                {
                    vals[i + 0 * ncols_m] = -vals_m[i];
                    vals[i + 1 * ncols_m] = vals_m[i] + dt * par[8] * vals_k[i];
                    col[i + 0 * ncols_m] = cols_m[i] + ind_var * nPoint +
                                           (ind_time - 1) * nPoint * state_num;
                    col[i + 1 * ncols_m] =
                        cols_m[i] + ind_var * nPoint + ind_time * nPoint * state_num;
                }
                MatSetValues(A, 1, &row, 2 * ncols_m, col, vals, INSERT_VALUES);

                break;
            }
        }
        MatRestoreRow(M, ind_point, NULL, NULL, &vals_m);
        MatRestoreRow(K, ind_point, NULL, NULL, &vals_k);
        MatRestoreRow(P[0], ind_point, NULL, NULL, &vals_px);
        MatRestoreRow(P[1], ind_point, NULL, NULL, &vals_py);

        PetscFree(col);
        PetscFree(vals);
    }
    // delete col, vals;
    // delete vals_m, vals_k, vals_px, vals_py;
}

void DiffConv2D::FormMatrixA21(Mat M, Mat K, Mat &A)
{
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void DiffConv2D::FormMatrixA22(Mat M, Mat K, Mat &A)
{
    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);

    // IS is_row, is_col;
    // Mat Msub;

    for (row = start; row < end; row++)
    {
        PetscInt ind_point, ind_var, ind_time;
        PetscReal scale;
        PetscInt *col;
        PetscReal *vals;
        PetscInt ncols_m;
        const PetscReal *vals_m;
        const PetscInt *cols_m;

        GetMatrixPosition(row, ctrl_num, ind_point, ind_var, ind_time);
        // MatGetRow(M,ind_point,&ncols_m, &cols_m,&vals_m);

        IS is_row, is_col;
        IS *is_row_array, *is_col_array;
        Mat *Msub;

        PetscMalloc1(1, &is_row_array);
        PetscMalloc1(1, &is_col_array);
        PetscMalloc1(1, &Msub);

        ISCreateStride(PETSC_COMM_SELF, 1, ind_point, 1, &is_row);
        ISCreateStride(PETSC_COMM_SELF, nPoint, 0, 1, &is_col);
        is_row_array[0] = is_row;
        is_col_array[0] = is_col;

        MatCreateSubMatrices(M, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Msub);
        MatGetRow(Msub[0], 0, &ncols_m, &cols_m, &vals_m);
        PetscMalloc1(ncols_m, &col);
        PetscMalloc1(ncols_m, &vals);

        // cout << "Ncols_k: " << ncols_k <<endl;
        // for(int i =0;i<ncols_k;i++)
        // 	cout << "Col: "<< cols_k[i] << " Value: " << vals_k[i] <<endl;
        // cout << "Ncols_m: " << ncols_m <<endl;
        // for(int i =0;i<ncols_m;i++)
        // 	cout << "Col: "<< cols_m[i] << " Value: " << vals_m[i] <<endl;

        // cout << "Output M_row, K_row Done!\n";

        switch (ind_var)
        {
        case 0:
            scale = beta1;
            break;
        case 1:
            scale = beta1;
            break;
        }

        switch (time_int)
        {
        case 0:
            // !steady state
            scale = scale;
            break;
        case 1:
            // !trapezoidal rule
            if (ind_time == 0 || ind_time == (nTstep - 1))
                scale = scale * 0.5 * dt;
            else
                scale = scale * dt;
            break;
        case 2:
            // !rectangle rule
            if (ind_time == (nTstep - 1))
                scale = 0.0;
            else
                scale = scale * dt;
            break;
        }

        for (int i = 0; i < ncols_m; i++)
        {
            vals[i] = vals_m[i] * scale;
            col[i] = cols_m[i] + ind_var * nPoint + ind_time * nPoint * ctrl_num;
        }

        // for(int i =0;i<ncols_m;i++)
        // 	cout << "Col: "<< col[i] << " Value: " << vals[i] <<endl;

        MatSetValues(A, 1, &row, ncols_m, col, vals, INSERT_VALUES);
        // MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
        MatRestoreRow(Msub[0], 0, NULL, NULL, &vals_m);

        PetscFree(is_row_array);
        PetscFree(is_col_array);
        PetscFree(Msub);
        PetscFree(col);
        PetscFree(vals);
    }
    // MatDestroy(&Msub);
}

void DiffConv2D::FormMatrixA23(Mat M, Mat K, Mat &A)
{
    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);

    for (row = start; row < end; row++)
    {
        PetscInt ind_point, ind_var, ind_time;
        PetscReal scale;
        PetscInt *col;
        PetscReal *vals;
        PetscInt ncols_m;
        const PetscReal *vals_m;
        const PetscInt *cols_m;

        GetMatrixPosition(row, ctrl_num, ind_point, ind_var, ind_time);
        // MatGetRow(M,ind_point,&ncols_m, &cols_m,&vals_m);

        IS is_row, is_col;
        IS *is_row_array, *is_col_array;
        Mat *Msub;

        PetscMalloc1(1, &is_row_array);
        PetscMalloc1(1, &is_col_array);
        PetscMalloc1(1, &Msub);

        ISCreateStride(PETSC_COMM_SELF, 1, ind_point, 1, &is_row);
        ISCreateStride(PETSC_COMM_SELF, nPoint, 0, 1, &is_col);
        is_row_array[0] = is_row;
        is_col_array[0] = is_col;

        MatCreateSubMatrices(M, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Msub);
        MatGetRow(Msub[0], 0, &ncols_m, &cols_m, &vals_m);

        // ISCreateStride(PETSC_COMM_WORLD, 1, ind_point, 0, &is_row);
        // ISCreateStride(PETSC_COMM_WORLD, nPoint, 0, 1, &is_col);
        // MatCreateSubMatrix(M, is_row, is_col, MAT_INITIAL_MATRIX, &Msub);
        // MatGetRow(Msub, 0, &ncols_m, &cols_m, &vals_m);

        PetscMalloc1(ncols_m, &col);
        PetscMalloc1(ncols_m, &vals);

        switch (time_int)
        {
        case 0:
            // !steady state
            scale = -1.0;
            break;
        case 1:
            // !trapezoidal rule
            scale = -dt;
            break;
        case 2:
            // !rectangle rule
            scale = -dt;
            break;
        }

        for (int i = 0; i < ncols_m; i++)
        {
            vals[i] = vals_m[i] * scale;
            col[i] = cols_m[i] + ind_var * nPoint + ind_time * nPoint * state_num;
        }

        MatSetValues(A, 1, &row, ncols_m, col, vals, INSERT_VALUES);

        // MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
        MatRestoreRow(Msub[0], 0, NULL, NULL, &vals_m);
        PetscFree(is_row_array);
        PetscFree(is_col_array);
        PetscFree(Msub);
        PetscFree(col);
        PetscFree(vals);
    }
}

void DiffConv2D::FormMatrixA31(Mat M, Mat K, Mat P[dim], Mat &A)
{
    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);

    IS is_row, is_col;
    Mat Msub, Ksub, Pxsub, Pysub, Convsub;

    for (row = start; row < end; row++)
    {
        PetscInt ind_point, ind_var, ind_time;
        PetscInt ncols_m, ncols_k, ncols_px, ncols_py, ncols_conv;
        const PetscReal *vals_m, *vals_k, *vals_px, *vals_py, *vals_conv;
        const PetscInt *cols_m, *cols_k, *cols_px, *cols_py, *cols_conv;

        PetscInt *col;
        PetscReal *vals;
        PetscReal scale;

        // cout << "Output M_row, K_row Done!\n";

        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
        // MatGetRow(M,ind_point,&ncols_m, &cols_m,&vals_m);
        // MatGetRow(K,ind_point,&ncols_k, &cols_k,&vals_k);
        // MatGetRow(P[0],ind_point,&ncols_px, &cols_px,&vals_px);
        // MatGetRow(P[1],ind_point,&ncols_py, &cols_py,&vals_py);
        // // ! Get Convection matrix
        // MatGetRow(Conv,ind_point,&ncols_conv, &cols_conv,&vals_conv);
        IS is_row, is_col;
        IS *is_row_array, *is_col_array;
        Mat *Msub, *Ksub, *Pxsub, *Pysub, *Convsub;

        PetscMalloc1(1, &is_row_array);
        PetscMalloc1(1, &is_col_array);
        PetscMalloc1(1, &Msub);
        PetscMalloc1(1, &Ksub);
        PetscMalloc1(1, &Pxsub);
        PetscMalloc1(1, &Pysub);
        PetscMalloc1(1, &Convsub);

        ISCreateStride(PETSC_COMM_SELF, 1, ind_point, 1, &is_row);
        ISCreateStride(PETSC_COMM_SELF, nPoint, 0, 1, &is_col);
        is_row_array[0] = is_row;
        is_col_array[0] = is_col;

        MatCreateSubMatrices(M, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Msub);
        MatCreateSubMatrices(K, 1, is_row_array, is_col_array, MAT_INITIAL_MATRIX,
                             &Ksub);
        MatCreateSubMatrices(P[0], 1, is_row_array, is_col_array,
                             MAT_INITIAL_MATRIX, &Pxsub);
        MatCreateSubMatrices(P[1], 1, is_row_array, is_col_array,
                             MAT_INITIAL_MATRIX, &Pysub);
        MatCreateSubMatrices(Conv, 1, is_row_array, is_col_array,
                             MAT_INITIAL_MATRIX, &Convsub);

        MatGetRow(Msub[0], 0, &ncols_m, &cols_m, &vals_m);
        MatGetRow(Ksub[0], 0, &ncols_k, &cols_k, &vals_k);
        MatGetRow(Pxsub[0], 0, &ncols_px, &cols_px, &vals_px);
        MatGetRow(Pysub[0], 0, &ncols_py, &cols_py, &vals_py);
        // ! Get Convection matrix
        MatGetRow(Convsub[0], 0, &ncols_conv, &cols_conv, &vals_conv);

        // cout << "Ncols_k: " << ncols_k <<endl;
        // for(int i =0;i<ncols_k;i++)
        // 	cout << "Col: "<< cols_k[i] << " Value: " << vals_k[i] <<endl;
        // cout << "Ncols_m: " << ncols_m <<endl;
        // for(int i =0;i<ncols_m;i++)
        // 	cout << "Col: "<< cols_m[i] << " Value: " << vals_m[i] <<endl;

        // cout << "ind_time: " << ind_time << " ind_point: " << ind_point<< "
        // ind_var: " << ind_var <<"\n";
        if (time_int != 0)
        {
            if (ind_time > 0)
            {
                switch (ind_var)
                {
                case 0:
                    PetscMalloc1(2 * ncols_m, &col);
                    PetscMalloc1(2 * ncols_m, &vals);
                    for (int i = 0; i < ncols_m; i++)
                    {
                        vals[i + 0 * ncols_m] = -vals_m[i];
                        vals[i + 1 * ncols_m] = vals_m[i] + dt * par[0] * vals_k[i];
                        col[i + 0 * ncols_m] =
                            cols_m[i] + 0 * nPoint + (ind_time - 1) * nPoint * state_num;
                        col[i + 1 * ncols_m] =
                            cols_m[i] + 0 * nPoint + ind_time * nPoint * state_num;
                    }
                    MatSetValues(A, 1, &row, ncols_m * 2, col, vals, INSERT_VALUES);

                    break;
                case 1:
                    PetscMalloc1(2 * ncols_m, &col);
                    PetscMalloc1(2 * ncols_m, &vals);
                    for (int i = 0; i < ncols_m; i++)
                    {
                        vals[i + 0 * ncols_m] = -vals_m[i];
                        vals[i + 1 * ncols_m] = vals_m[i] + dt * par[0] * vals_k[i];
                        col[i + 0 * ncols_m] =
                            cols_m[i] + 1 * nPoint + (ind_time - 1) * nPoint * state_num;
                        col[i + 1 * ncols_m] =
                            cols_m[i] + 1 * nPoint + ind_time * nPoint * state_num;
                    }
                    MatSetValues(A, 1, &row, (2 * ncols_m), col, vals, INSERT_VALUES);
                    break;
                }
            }
            else
            {
                switch (ind_var)
                {
                case 0:
                    PetscMalloc1(1 * ncols_m, &col);
                    PetscMalloc1(1 * ncols_m, &vals);
                    for (int i = 0; i < ncols_m; i++)
                    {
                        vals[i] = vals_m[i] + dt * par[0] * vals_k[i];
                        col[i] = cols_m[i] + 0 * nPoint + ind_time * nPoint * state_num;
                    }
                    MatSetValues(A, 1, &row, ncols_m, col, vals, INSERT_VALUES);

                    break;
                case 1:
                    PetscMalloc1(1 * ncols_m, &col);
                    PetscMalloc1(1 * ncols_m, &vals);
                    for (int i = 0; i < ncols_m; i++)
                    {
                        vals[i] = vals_m[i] + dt * par[0] * vals_k[i];
                        col[i] = cols_m[i] + 1 * nPoint + ind_time * nPoint * state_num;
                    }
                    MatSetValues(A, 1, &row, ncols_m, col, vals, INSERT_VALUES);
                    break;
                }
            }
        }
        else
        {
            PetscMalloc1(1 * ncols_k, &col);
            PetscMalloc1(1 * ncols_k, &vals);
            for (int i = 0; i < ncols_k; i++)
            {
                // ! Pure Diffusion
                // vals[i] = par[0] * vals_k[i];
                // ! Convection
                vals[i] = par[0] * vals_k[i] + vals_conv[i];

                col[i] = cols_k[i];
            }
            MatSetValues(A, 1, &row, ncols_k, col, vals, INSERT_VALUES);
        }

        // MatRestoreRow(M,ind_point,NULL,NULL,&vals_m);
        // MatRestoreRow(K,ind_point,NULL,NULL,&vals_k);
        // MatRestoreRow(P[0],ind_point,NULL,NULL,&vals_px);
        // MatRestoreRow(P[1],ind_point,NULL,NULL,&vals_py);
        // MatRestoreRow(Conv,ind_point,NULL,NULL,&vals_conv);
        MatRestoreRow(Msub[0], 0, NULL, NULL, &vals_m);
        MatRestoreRow(Ksub[0], 0, NULL, NULL, &vals_k);
        MatRestoreRow(Pxsub[0], 0, NULL, NULL, &vals_px);
        MatRestoreRow(Pysub[0], 0, NULL, NULL, &vals_py);
        MatRestoreRow(Convsub[0], 0, NULL, NULL, &vals_conv);

        PetscFree(is_row_array);
        PetscFree(is_col_array);
        PetscFree(Msub);
        PetscFree(Ksub);
        PetscFree(Pxsub);
        PetscFree(Pysub);
        PetscFree(Convsub);
        PetscFree(col);
        PetscFree(vals);
    }
}

// ! Not used, compute A23 and transpose to get A32
void DiffConv2D::FormMatrixA32(Mat M, Mat K, Mat &A)
{
    PetscInt ind_point, ind_var, ind_time;
    PetscInt row, start, end;
    PetscInt *col = new PetscInt[nPoint];
    PetscReal *vals = new PetscReal[nPoint];
    const PetscReal *vals_m;
    PetscReal scale;
    MatGetOwnershipRange(A, &start, &end);

    for (row = start; row < end; row++)
    {
        GetMatrixPosition(row, state_num, ind_point, ind_var, ind_time);
        MatGetRow(M, ind_point, NULL, NULL, &vals_m);
        scale = dt;
        for (int i = 0; i < nPoint; i++)
        {
            vals[i] = vals_m[i] * scale;
            col[i] = i + ind_var * nPoint + ind_time * nPoint * ctrl_num;
        }
        MatSetValues(A, 1, &row, nPoint, col, vals, INSERT_VALUES);
        MatRestoreRow(M, ind_point, NULL, NULL, &vals_m);
    }

    delete col;
    delete vals;
}

void DiffConv2D::FormMatrixA33(Mat M, Mat K, Mat &A)
{
    // explicitly set diagonal to be zero
    PetscInt row, start, end;
    MatGetOwnershipRange(A, &start, &end);
    for (row = start; row < end; row++)
        MatSetValue(A, row, row, 0.0, INSERT_VALUES);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}

void DiffConv2D::ApplyBoundaryCondition(const UserSetting2D *ctx, int flag)
{
    PetscInt *bc_global_ids;
    PetscScalar *bc_vals;

    // Vec Rhs_BC;
    // VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, (state_num * 2 + ctrl_num) *
    // nPoint * nTstep, &Rhs_BC); VecSet(Rhs_BC, 0.0);

    PetscMalloc1(nTstep * n_bcval * 2, &bc_global_ids);
    PetscMalloc1(nTstep * n_bcval * 2, &bc_vals);
    int count = 0;
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {
            if (ctx->bc_flag[i] == -1)
            {
                bc_global_ids[count + 0] = i + 0 * nPoint + j * state_num * nPoint;
                // bc_vals[count + 0] = ctx->val_bc[0][i]-State[0][i + j * nPoint]; //!
                // bc for nonlinear iteration
                bc_vals[count + 0] = ctx->val_bc[0][i]; //! bc for linear
                count += 1;
                continue;
            }
        }
    }

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {

            if (ctx->bc_flag[i] == -1)
            {
                bc_global_ids[count + 0] = i + 0 * nPoint + j * state_num * nPoint +
                                           (state_num + ctrl_num) * nPoint * nTstep;
                // bc_global_ids[count + 1] = i + 1 * nPoint + j * state_num * nPoint +
                // (state_num+ctrl_num) *nPoint *nTstep; bc_vals[count + 0] = 0.0 -
                // Lambda[0][i + j * nPoint];
                bc_vals[count + 0] = 0.0; //! bc for linear
                // bc_vals[count + 1] = 0.0 - Lambda[1][i + j * nPoint];
                count += 1;
                continue;
            }
        }
    }

    // string fname = work_dir + "debug/bc_global_ids.txt";
    // ofstream fout;
    // fout.open(fname.c_str());

    // if (fout.is_open())
    // {
    //     for (int i = 0; i < nTstep * n_bcval; i++)
    //         fout << bc_global_ids[i] << endl;
    // }
    // else
    // {
    //     cout << "Cannot open " << fname << "!\n";
    // }

    if (flag == 0)
    {
        MatZeroRowsColumns(TanMat, nTstep * n_bcval * 2, bc_global_ids, 1.0, 0, 0);
        PetscPrintf(PETSC_COMM_WORLD, "Apply BC in Tangent matrix Done!\n");
    }
    MPI_Barrier(comm);

    VecSetValues(ResVec, nTstep * n_bcval * 2, bc_global_ids, bc_vals,
                 INSERT_VALUES); // ! Modify to test boundary condition
    VecAssemblyBegin(ResVec);
    VecAssemblyEnd(ResVec);
    PetscPrintf(PETSC_COMM_WORLD, "Apply BC in Residual vector Done!\n");

    PetscFree(bc_vals);
    PetscFree(bc_global_ids);
    delete bc_vals;
    delete bc_global_ids;
}

void DiffConv2D::ApplyBoundaryConditionNest(const UserSetting2D *ctx,
                                            int flag)
{
    PetscInt *bc_global_ids;
    PetscScalar *bc_vals;

    // Vec Rhs_BC;
    // VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, (state_num * 2 + ctrl_num) *
    // nPoint * nTstep, &Rhs_BC); VecSet(Rhs_BC, 0.0);

    PetscMalloc1(nTstep * n_bcval * 2, &bc_global_ids);
    PetscMalloc1(nTstep * n_bcval * 2, &bc_vals);
    int count = 0;
    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {
            if (ctx->bc_flag[i] == -1)
            {
                bc_global_ids[count + 0] = i + 0 * nPoint + j * state_num * nPoint;
                // bc_vals[count + 0] = ctx->val_bc[0][i]-State[0][i + j * nPoint]; //!
                // bc for nonlinear iteration
                bc_vals[count + 0] = ctx->val_bc[0][i]; //! bc for linear
                count += 1;
                continue;
            }
        }
    }

    for (int j = 0; j < nTstep; j++)
    {
        for (int i = 0; i < ctx->bc_flag.size(); i++)
        {

            if (ctx->bc_flag[i] == -1)
            {
                bc_global_ids[count + 0] = i + 0 * nPoint + j * state_num * nPoint +
                                           (state_num + ctrl_num) * nPoint * nTstep;
                // bc_global_ids[count + 1] = i + 1 * nPoint + j * state_num * nPoint +
                // (state_num+ctrl_num) *nPoint *nTstep; bc_vals[count + 0] = 0.0 -
                // Lambda[0][i + j * nPoint];
                bc_vals[count + 0] = 0.0; //! bc for linear
                // bc_vals[count + 1] = 0.0 - Lambda[1][i + j * nPoint];
                count += 1;
                continue;
            }
        }
    }

    // string fname = work_dir + "debug/bc_global_ids.txt";
    // ofstream fout;
    // fout.open(fname.c_str());

    // if (fout.is_open())
    // {
    //     for (int i = 0; i < nTstep * n_bcval; i++)
    //         fout << bc_global_ids[i] << endl;
    // }
    // else
    // {
    //     cout << "Cannot open " << fname << "!\n";
    // }

    if (flag == 0)
    {
        MatZeroRowsColumns(TanMat, nTstep * n_bcval * 2, bc_global_ids, 1.0, 0, 0);
        PetscPrintf(PETSC_COMM_WORLD, "Apply BC in Tangent matrix Done!\n");
    }

    VecSetValues(ResVec, nTstep * n_bcval * 2, bc_global_ids, bc_vals,
                 INSERT_VALUES); // ! Modify to test boundary condition
    VecAssemblyBegin(ResVec);
    VecAssemblyEnd(ResVec);
    PetscPrintf(PETSC_COMM_WORLD, "Apply BC in Residual vector Done!\n");

    PetscFree(bc_vals);
    PetscFree(bc_global_ids);
    delete bc_vals;
    delete bc_global_ids;
}

void DiffConv2D::TangentMatSetup()
{

    /// A11
    FormMatrixA11(M, K, Asubmat[0]); // need debug
    MatAssemblyBegin(Asubmat[0], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[0], MAT_FINAL_ASSEMBLY);

    DebugVisualizeMat(Asubmat[0], work_dir + "debug/matA11.m");

    // A12 & A21 // No problem
    FormMatrixA12(M, K, Asubmat[1]); // need debug
    MatAssemblyBegin(Asubmat[1], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[1], MAT_FINAL_ASSEMBLY);
    MatTranspose(Asubmat[1], MAT_INITIAL_MATRIX, &Asubmat[3]); // correct

    DebugVisualizeMat(Asubmat[1], work_dir + "debug/matA12.m");
    DebugVisualizeMat(Asubmat[3], work_dir + "debug/matA21.m");

    // A13 & A31 // No problem
    FormMatrixA31(M, K, P, Asubmat[6]); // need debug
    MatAssemblyBegin(Asubmat[6], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[6], MAT_FINAL_ASSEMBLY);
    //! Stablization mat
    MatAXPY(Asubmat[6], 1.0, Stable, DIFFERENT_NONZERO_PATTERN);
    // PetscObjectSetName((PetscObject) Stable, "MatStable");
    // DebugVisualizeMat(Stable, work_dir + "debug/Stable.m");

    MatTranspose(Asubmat[6], MAT_INITIAL_MATRIX, &Asubmat[2]); // correct

    DebugVisualizeMat(Asubmat[2], work_dir + "debug/matA13.m");
    DebugVisualizeMat(Asubmat[6], work_dir + "debug/matA31.m");

    // A22  // No problem
    FormMatrixA22(M, K, Asubmat[4]);
    MatAssemblyBegin(Asubmat[4], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[4], MAT_FINAL_ASSEMBLY);

    DebugVisualizeMat(Asubmat[4], work_dir + "debug/matA22.m");

    // A23 & A32 // No problem
    FormMatrixA23(M, K, Asubmat[5]);
    MatAssemblyBegin(Asubmat[5], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[5], MAT_FINAL_ASSEMBLY);
    MatTranspose(Asubmat[5], MAT_INITIAL_MATRIX, &Asubmat[7]);

    DebugVisualizeMat(Asubmat[5], work_dir + "debug/matA23.m");
    DebugVisualizeMat(Asubmat[7], work_dir + "debug/matA32.m");

    // A33 // No problem
    FormMatrixA33(M, K, Asubmat[8]); // need debug
    MatAssemblyBegin(Asubmat[8], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Asubmat[8], MAT_FINAL_ASSEMBLY);

    DebugVisualizeMat(Asubmat[8], work_dir + "debug/matA33.m");
    PetscBarrier((PetscObject)Asubmat);

    PetscPrintf(PETSC_COMM_WORLD, "Compute Submatrices Done!\n");

    MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, Asubmat, &TanMat_tmp);
    // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, work_dir +
    // "debug/TanMat_tmp.m",&viewer); ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(TanMat_tmp,viewer);

    PetscPrintf(PETSC_COMM_WORLD, "Create Mat Nest Done!\n");

    MatConvert(TanMat_tmp, MATMPIAIJ, MAT_INITIAL_MATRIX, &TanMat);
    PetscObjectSetName((PetscObject)TanMat, "TanMat");

    PetscPrintf(PETSC_COMM_WORLD, "Tangent Mat Convert Done!\n");
    MatSetUp(TanMat);

    // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, work_dir +
    // "debug/TanMat.m",&viewer); ierr = PetscViewerPushFormat(viewer,format);
    // ierr = MatView(TanMat,viewer);

    DebugVisualizeMat(TanMat, work_dir + "debug/TanMat_beforeBC.m");

    PetscPrintf(PETSC_COMM_WORLD, "finish Tangent Mat\n");
}

void DiffConv2D::TangentMatSetupAsubmat()
{

    /// A11

    // DebugVisualizeMat(Asubmat[0], work_dir + "debug/matA11.m");
    // // getchar();
    // cout << "A11 View done!\n";

    // A12 & A21 // No problem
    FormMatrixA12(M, K, Asubmat[1]); // need debug
    // MatAssemblyBegin(Asubmat[1], MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(Asubmat[1], MAT_FINAL_ASSEMBLY);
    MatTranspose(Asubmat[1], MAT_INITIAL_MATRIX, &Asubmat[3]); // correct

    // DebugVisualizeMat(Asubmat[1], work_dir + "debug/matA12.m");
    // // getchar();
    // cout << "A12 View done!\n";
    // DebugVisualizeMat(Asubmat[3], work_dir + "debug/matA21.m");
    // // getchar();
    // cout << "A21 View done!\n";

    // A13 & A31 // No problem
    MatTranspose(Asubmat[6], MAT_INITIAL_MATRIX, &Asubmat[2]); // correct

    // DebugVisualizeMat(Asubmat[2], work_dir + "debug/matA13.m");
    // // getchar();
    // cout << "A13 View done!\n";
    // DebugVisualizeMat(Asubmat[6], work_dir + "debug/matA31.m");
    // // getchar();
    // cout << "A31 View done!\n";

    // A22  // No problem
    // DebugVisualizeMat(Asubmat[4], work_dir + "debug/matA22.m");
    // // getchar();
    // cout << "A22 View done!\n";

    // A23 & A32 // No problem
    MatTranspose(Asubmat[5], MAT_INITIAL_MATRIX, &Asubmat[7]);

    // DebugVisualizeMat(Asubmat[5], work_dir + "debug/matA23.m");
    // // getchar();
    // cout << "A23 View done!\n";
    // DebugVisualizeMat(Asubmat[7], work_dir + "debug/matA32.m");
    // // getchar();
    // cout << "A32 View done!\n";
    // A33 // No problem
    FormMatrixA33(M, K, Asubmat[8]);

    // DebugVisualizeMat(Asubmat[8], work_dir + "debug/matA33.m");
    // // getchar();
    // cout << "A33 View done!\n";

    PetscPrintf(PETSC_COMM_WORLD, "Compute Submatrices Done!\n");

    MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, Asubmat, &TanMat_tmp);
    // DebugVisualizeMat(TanMat_tmp, work_dir + "debug/TanMat_tmp.m");

    PetscPrintf(PETSC_COMM_WORLD, "Create Mat Nest Done!\n");

    MatNestGetISs(TanMat_tmp, isg, NULL);
    DebugVisualizeIS(isg[0], work_dir + "debug/isg0_nest.m");
    DebugVisualizeIS(isg[1], work_dir + "debug/isg1_nest.m");
    DebugVisualizeIS(isg[2], work_dir + "debug/isg2_nest.m");

    // // * MATNEST convert to MATAIJ
    // Mat TanMat;
    // MatConvert(TanMat_tmp, MATMPIAIJ, MAT_INITIAL_MATRIX, &TanMat);
    // PetscObjectSetName((PetscObject)TanMat, "TanMat");

    // PetscPrintf(PETSC_COMM_WORLD, "Tangent Mat Convert Done!\n");

    // DebugVisualizeMat(TanMatConvert, work_dir +
    //                                      "debug/TanMatConvert_beforeBC.m");

    PetscPrintf(PETSC_COMM_WORLD, "finish Tangent Mat\n");
}

void DiffConv2D::FormDesireVec(const UserSetting2D *user)
{
    // Build linear system in each process
    int e;
    PetscPrintf(PETSC_COMM_WORLD, "Computing Residual Vector.\n");
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        int nen, A;

        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;

        vector<double> qtmp;
        double tmp_desire[state_num];

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0.0);
        dNdx.clear();
        dNdx.resize(nen, {0.0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0.0}});
        qtmp.clear();
        qtmp.resize(nen * state_num * nTstep, 0.0);

        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                double X(0.0), Y(0.0), Z(0.0), T(0.0);
                double pt[2];
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;
                bzmesh_process[e].Para2Phys(Gpt[i], Gpt[j], pt);
                X = pt[0];
                Y = pt[1];
                for (int tt = 0; tt < nTstep; tt++)
                {
                    T = tt * dt;
                    for (int k = 0; k < nen; k++)
                    {
                        user->DesireStateFunction(X, Y, Z, T, tmp_desire);
                        for (int s = 0; s < state_num; s++)
                            qtmp[tt * state_num * nen + s * nen + k] +=
                                tmp_desire[s] * Nx[k] * detJ;
                    }
                }
            }
        }

        ResnlAssembly(qtmp, bzmesh_process[e].IEN, Res_Desire);
    }
    cout << "Process " << comRank << " :complete compute Desire vector!\n";
    VecAssemblyBegin(Res_Desire);
}

void DiffConv2D::FormDesireVec_Vertex(const UserSetting2D *user)
{
    // Build linear system in each process
    int e;
    PetscPrintf(PETSC_COMM_WORLD, "Computing Residual Vector.\n");
    cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "
         << bzmesh_process.size() << " elements.\n";
    for (e = 0; e < bzmesh_process.size(); e++)
    {
        int nen, A;

        double dudx[dim][dim];
        double detJ;
        vector<double> Nx;
        vector<array<double, 2>> dNdx;
        vector<array<array<double, 2>, 2>> dN2dx2;

        vector<double> qtmp;
        double tmp_desire[state_num];

        nen = bzmesh_process[e].IEN.size();

        Nx.clear();
        Nx.resize(nen, 0.0);
        dNdx.clear();
        dNdx.resize(nen, {0.0});
        dN2dx2.clear();
        dN2dx2.resize(nen, {{0.0}});
        qtmp.clear();
        qtmp.resize(nen * state_num * nTstep, 0.0);

        for (int i = 0; i < Gpt.size(); i++)
        {
            for (int j = 0; j < Gpt.size(); j++)
            {
                double X(0.0), Y(0.0), Z(0.0), T(0.0);
                BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts,
                              bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
                detJ = wght[i] * wght[j] * detJ;

                for (int tt = 0; tt < nTstep; tt++)
                {
                    T = tt * dt;
                    for (int k = 0; k < nen; k++)
                    {

                        for (int kk = 0; kk < nen; kk++)
                        {
                            for (int s = 0; s < state_num; s++)
                                qtmp[tt * state_num * nen + s * nen + k] +=
                                    user->val_desire[s][bzmesh_process[e].IEN[kk]] * Nx[kk] * Nx[k] * detJ;
                        }
                    }
                }
            }
        }

        ResnlAssembly(qtmp, bzmesh_process[e].IEN, Res_Desire);
    }
    cout << "Process " << comRank << " :complete compute Desire vector!\n";
    VecAssemblyBegin(Res_Desire);
}

void DiffConv2D::FormResVecb1(Vec x)
{
    Vec y_k, l_k, u_k;
    Vec vec_tmp0, vec_tmp1, vec_tmp2;

    // MatNestGetISs(TanMat_tmp, isg, NULL);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                0, 1, &isg[0]);
    // ISCreateStride(PETSC_COMM_WORLD, ctrl_num * nPoint * nTstep,
    //                state_num * nPoint * nTstep, 1, &isg[1]);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                (state_num + ctrl_num) * nPoint * nTstep, 1,
    //                &isg[2]);

    VecGetSubVector(x, isg[0], &y_k);
    VecGetSubVector(x, isg[1], &u_k);
    VecGetSubVector(x, isg[2], &l_k);

    // DebugVisualizeIS(isg[0], work_dir + "debug/isg0.m");
    // DebugVisualizeIS(isg[1], work_dir + "debug/isg1.m");
    // DebugVisualizeIS(isg[2], work_dir + "debug/isg2.m");

    VecDuplicate(y_k, &bsub[0]);
    VecDuplicate(y_k, &vec_tmp0);
    VecDuplicate(y_k, &vec_tmp1);
    VecDuplicate(y_k, &vec_tmp2);

    // ! Newton's method for nonlinear
    // VecWAXPY(vec_tmp0, -1.0, y_k, Y_d);		 //-(y_k - y_d)
    // MatMult(Asubmat[0], vec_tmp0, vec_tmp1); //vec_tmp1 = -A11 * (y_k-y_d)

    // MatMult(Asubmat[2], l_k, vec_tmp2); // vec_tmp2 = A13 * l_k
    // VecWAXPY(bsub[0], -1.0, vec_tmp2, vec_tmp1); //vec_tmp1 - vec_tmp2 = -A11 *
    // (y_k-y_d) - A13 * l_k

    // ! Heat equation
    // MatMult(Asubmat[0], Y_d, bsub[0]); // * bsub[0] = A11 * y_d
    VecCopy(Res_Desire, bsub[0]); // * bsub[0] = Res_Desire

    // ! Apply bc to vector
    Vec vec_bc;
    VecDuplicate(Y_bc, &vec_bc);
    MatMult(Asubmat[0], Y_bc, vec_bc);
    VecAXPY(bsub[0], -1, vec_bc);

    PetscObjectSetName((PetscObject)bsub[0], "bsub1");

    VecRestoreSubVector(x, isg[0], &y_k);
    VecRestoreSubVector(x, isg[1], &u_k);
    VecRestoreSubVector(x, isg[2], &l_k);
}

void DiffConv2D::FormResVecb2(Vec x)
{
    Vec y_k, l_k, u_k;
    Vec vec_tmp0, vec_tmp1, vec_tmp2;

    // MatNestGetISs(TanMat_tmp, isg, NULL);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                0, 1, &isg[0]);
    // ISCreateStride(PETSC_COMM_WORLD, ctrl_num * nPoint * nTstep,
    //                state_num * nPoint * nTstep, 1, &isg[1]);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                (state_num + ctrl_num) * nPoint * nTstep, 1,
    //                &isg[2]);

    VecGetSubVector(x, isg[0], &y_k);
    VecGetSubVector(x, isg[1], &u_k);
    VecGetSubVector(x, isg[2], &l_k);

    VecDuplicate(u_k, &bsub[1]);
    VecDuplicate(u_k, &vec_tmp0);
    VecDuplicate(u_k, &vec_tmp1);
    VecDuplicate(u_k, &vec_tmp2);

    // ! Newton's method for nonlinear
    // MatMult(Asubmat[4], u_k, vec_tmp1); //A22 * u_k
    // VecScale(vec_tmp1, -1.0); // vec_tmp1 = -A22 * u_k

    // MatMult(Asubmat[5], l_k, vec_tmp2); // vec_tmp2 = A13 * l_k
    // VecWAXPY(bsub[1], -1.0, vec_tmp2, vec_tmp1); //vec_tmp1 - vec_tmp2 = -A22 *
    // u_k - A23 * l_k

    // ! Heat equation
    VecSet(bsub[1], 0.0);
    PetscObjectSetName((PetscObject)bsub[1], "bsub2");

    VecRestoreSubVector(x, isg[0], &y_k);
    VecRestoreSubVector(x, isg[1], &u_k);
    VecRestoreSubVector(x, isg[2], &l_k);
}

void DiffConv2D::FormResVecb3(Vec x)
{

    Vec y_k, l_k, u_k;
    Vec vec_tmp0, vec_tmp1, vec_tmp2, vec_tmp3;

    // MatNestGetISs(TanMat_tmp, isg, NULL);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                0, 1, &isg[0]);
    // ISCreateStride(PETSC_COMM_WORLD, ctrl_num * nPoint * nTstep,
    //                state_num * nPoint * nTstep, 1, &isg[1]);
    // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
    //                (state_num + ctrl_num) * nPoint * nTstep, 1,
    //                &isg[2]);

    VecGetSubVector(x, isg[0], &y_k);
    VecGetSubVector(x, isg[1], &u_k);
    VecGetSubVector(x, isg[2], &l_k);

    VecDuplicate(l_k, &vec_tmp0);
    VecDuplicate(l_k, &vec_tmp1);
    VecDuplicate(l_k, &vec_tmp2);
    VecDuplicate(l_k, &vec_tmp3);
    VecDuplicate(l_k, &bsub[2]);

    // ! Newton's method for nonlinear
    // MatMult(Asubmat[6], y_k, vec_tmp1); //A31 * y_k
    // VecScale(vec_tmp1, -1.0); // vec_tmp1 = -A31 * y_k

    // MatMult(Asubmat[7], u_k, vec_tmp2); //A32 * u_k
    // VecScale(vec_tmp2, -1.0); // vec_tmp2 = -A32 * u_k

    // VecCopy(Res_nl, vec_tmp3);
    // // VecSet(vec_tmp3, 0.0);
    // VecWAXPY(vec_tmp0, 1.0, vec_tmp1, vec_tmp2); //vec_tmp1 + vec_tmp2 = -A31 *
    // y_k -A32 * u_k VecWAXPY(bsub[2], 1.0, vec_tmp3, vec_tmp0); //vec_tmp0 +
    // vec_tmp3 = d(y_k) -A31 * y_k -A32 * u_k

    // ! Heat equation
    VecSet(bsub[2], 0.0);
    // VecCopy(Res_nl, bsub[2]);

    // ! Apply bc to vector
    Vec vec_bc;
    VecDuplicate(Y_bc, &vec_bc);
    MatMult(Asubmat[6], Y_bc, vec_bc);
    VecAXPY(bsub[2], -1, vec_bc);

    PetscObjectSetName((PetscObject)bsub[2], "bsub3");

    VecRestoreSubVector(x, isg[0], &y_k);
    VecRestoreSubVector(x, isg[1], &u_k);
    VecRestoreSubVector(x, isg[2], &l_k);
}

void DiffConv2D::ResidualVecSetup(Vec x)
{
    FormResVecb1(x);
    DebugVisualizeVec(bsub[0], work_dir + "debug/VecB1.m");

    FormResVecb2(x);
    DebugVisualizeVec(bsub[1], work_dir + "debug/VecB2.m");

    FormResVecb3(x);
    DebugVisualizeVec(bsub[2], work_dir + "debug/VecB3.m");

    PetscInt low, high, ldim, iglobal, i, k;
    PetscReal val, *array;
    PetscInt shift[3] = {0, state_num * nTstep * nPoint,
                         (state_num + ctrl_num) * nTstep * nPoint};

    for (k = 0; k < 3; k++)
    {
        VecGetOwnershipRange(bsub[k], &low, &high);
        VecGetLocalSize(bsub[k], &ldim);
        VecGetArray(bsub[k], &array);
        // cout << "Rank: " << comRank <<" k: " << k << " low " << low << " high "
        // << high << "\n";
        for (i = 0; i < ldim; i++)
        {
            val = array[i];
            iglobal = i + low + shift[k];
            VecSetValues(ResVec, 1, &iglobal, &val, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(ResVec);
    VecAssemblyEnd(ResVec);

    PetscObjectSetName((PetscObject)ResVec, "ResVec");

    DebugVisualizeVec(ResVec, work_dir + "debug/ResVec_beforeBC.m");
    PetscPrintf(PETSC_COMM_WORLD, "finish Residual Vector\n");
}

void DiffConv2D::Run(const UserSetting2D *ctx)
{

    int n_iterator(1);
    double tol(1e-7);
    int l(0);
    time_t t0, t1;
    vector<double> Y_delta(state_num * nPoint * nTstep);
    vector<double> U_delta(ctrl_num * nPoint * nTstep);
    vector<double> L_delta(state_num * nPoint * nTstep);

    InitializeProblem(ctx);

    ReadBezierElementProcess(ctx->work_dir);

    PetscPrintf(PETSC_COMM_WORLD, "Building Unit Matrix...\n");
    t0 = time(NULL);

    // ? Assemble Submatrix directly, needs debugging
    bool flag_glb = true;

    // if(!flag_glb)
    {
        t0 = time(NULL);
        BuildLinearSystemProcessAsubmat(ctx->pts, ctx->val_bc, ctx->val_ini, false);
        MatAssemblyEnd(Asubmat[0], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Asubmat[4], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Asubmat[5], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Asubmat[6], MAT_FINAL_ASSEMBLY);
        TangentMatSetupAsubmat();
        t1 = time(NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Done Sub Matrix Assembly with time: %d \n",
                    t1 - t0);
    }
    // else

    {
        t0 = time(NULL);
        BuildLinearSystemProcessAsubmat(ctx->pts, ctx->val_bc, ctx->val_ini,
                                        flag_glb);
        MatAssemblyEnd(TanMat, MAT_FINAL_ASSEMBLY);
        t1 = time(NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n",
                    t1 - t0);
    }
    // DebugVisualizeMat(TanMat, work_dir + "debug/TanMat2_BeforeBC.m");

    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    if (!ctx->GetReadDesire())
        FormDesireVec(ctx);
    else
        FormDesireVec_Vertex(ctx);

    VecAssemblyEnd(Res_Desire);
    DebugVisualizeVec(Res_Desire, work_dir + "debug/DesireVec.m");

    t1 = time(NULL);
    PetscPrintf(PETSC_COMM_WORLD, "Done Unit Matrix Assembly with time: %d...\n",
                t1 - t0);
    MPI_Barrier(comm);
    // PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");
    // t0 = time(NULL);
    //  TangentMatSetup();
    // t1 = time(NULL);
    // PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 -
    // t0);

    // Old KSP code
    for (l = 0; l < n_iterator; l++)
    {
        stringstream ss;
        ss << l;
        PetscPrintf(PETSC_COMM_WORLD, "Iteration Step: %d\n", l);

        PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");
        // //if(!flag_glb)
        {
            t0 = time(NULL);
            // TangentMatSetup();
            // TangentMatSetupAsubmat(); // ? Assemble Submatrix directly, needs
            // debugging
            t1 = time(NULL);
            PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n",
                        t1 - t0);
        }

        MPI_Barrier(comm);
        VecSet(Res_nl, 0.0);
        BuildResVectorProcess(ctx->pts, ctx->val_bc, ctx->val_ini, ctx->bc_flag);
        VecAssemblyEnd(Res_nl);

        // DebugVisualizeMat(M,  work_dir + "debug/MMat.m");
        // DebugVisualizeMat(K,  work_dir + "debug/KMat.m");
        // DebugVisualizeMat(P[0],  work_dir + "debug/PxMat.m");
        // DebugVisualizeMat(P[1],  work_dir + "debug/PyMat.m");
        // DebugVisualizeMat(Conv, work_dir + "debug/ConvMat.m");

        ResidualVecSetup(X);
        MPI_Barrier(comm);

        //! Apply BC:Modify for stabilization method
        ApplyBoundaryCondition(ctx, l);
        // ApplyBoundaryCondition(ctx, 0);
        // DebugSubMat();
        // DebugVisualizeMat(PCMat,  work_dir + "debug/PCMat.m");
        //  DebugVisualizeMat(TanMat,  work_dir + "debug/TanMat_afterBC.m");
        //  DebugVisualizeVec(ResVec,  work_dir + "debug/ResVec_afterBC.m");

        // t1 = time(NULL);
        // PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1
        // - t0);
        MPI_Barrier(comm);
        PetscPrintf(PETSC_COMM_WORLD, "Solving...\n");

        t0 = time(NULL);

        /*Petsc KSP solver setting*/
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        // KSPSetOperators(ksp, TanMat, TanMat);
        KSPSetOperators(ksp, TanMat, TanMat);

        KSPSetFromOptions(ksp);
        KSPSolve(ksp, ResVec, dX);
        PetscInt its;
        KSPGetIterationNumber(ksp, &its);
        PetscPrintf(PETSC_COMM_WORLD, "iterations %d\n", its);
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
        t1 = time(NULL);
        PetscPrintf(PETSC_COMM_WORLD, "Done Solving with time: %d \n", t1 - t0);
        // KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
        PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");

        PetscPrintf(PETSC_COMM_WORLD, "Updating X...\n");
        PetscObjectSetName((PetscObject)dX, "dX");
        DebugVisualizeVec(X, work_dir + "debug/X" + ss.str() + ".m");
        DebugVisualizeVec(dX, work_dir + "debug/dX" + ss.str() + ".m");

        CollectAndUpdateResult(dX, X);

        VecAXPY(X, 1.0, dX);

        PetscPrintf(PETSC_COMM_WORLD, "Visualizing...\n");
        for (int i = 0; i < nTstep; i++)
        {
            // VisualizeVTK_ControlMesh_Burger(ctx->pts, ctx->mesh, i, l, work_dir +
            // "test1/");
            VisualizeVTK_ControlMesh_Heat(ctx->pts, ctx->mesh, i, l, ctx->output_dir);
            VisualizeVTK_PhysicalDomain(i, l, ctx->output_dir);
        }

        PetscReal ResVal;
        PetscReal ResVal_y, ResVal_u, ResVal_l;
        Vec dy, dl, du;

        // MatNestGetISs(TanMat_tmp, isg, NULL);
        // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
        //                0, 1, &isg[0]);
        // ISCreateStride(PETSC_COMM_WORLD, ctrl_num * nPoint * nTstep,
        //                state_num * nPoint * nTstep, 1, &isg[1]);
        // ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
        //                (state_num + ctrl_num) * nPoint * nTstep, 1,
        //                &isg[2]);

        VecGetSubVector(dX, isg[0], &dy);
        VecGetSubVector(dX, isg[1], &du);
        VecGetSubVector(dX, isg[2], &dl);
        VecNorm(dy, NORM_2, &ResVal_y);
        VecNorm(du, NORM_2, &ResVal_u);
        VecNorm(dl, NORM_2, &ResVal_l);

        VecNorm(ResVec, NORM_2, &ResVal);
        PetscPrintf(PETSC_COMM_WORLD, "Residual: %lf \n", ResVal);
        PetscPrintf(PETSC_COMM_WORLD, "|dy| = %f;|du| = %f; |dl|=%f \n", ResVal_y,
                    ResVal_u, ResVal_l);
        PetscPrintf(PETSC_COMM_WORLD, "|dy|+|du|+|dl| %lf \n",
                    ResVal_y + ResVal_u + ResVal_l);

        PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");

        PetscReal ObjFunVal, ObjFunVal_y, ObjFunVal_u;
        Vec y_new, u_new;
        Vec y_tmp, u_tmp;
        VecGetSubVector(X, isg[0], &y_new);
        VecGetSubVector(X, isg[1], &u_new);
        VecDuplicate(y_new, &y_tmp);
        VecDuplicate(u_new, &u_tmp);

        MatMult(Asubmat[0], y_new, y_tmp);
        VecDot(y_new, y_tmp, &ObjFunVal_y);
        PetscPrintf(PETSC_COMM_WORLD, "Object function y value: %lf \n",
                    ObjFunVal_y);

        MatMult(Asubmat[4], u_new, u_tmp);
        VecDot(u_new, u_tmp, &ObjFunVal_u);
        PetscPrintf(PETSC_COMM_WORLD, "Object function u value: %lf \n",
                    ObjFunVal_u);

        ObjFunVal = 0.5 * (ObjFunVal_y + ObjFunVal_u);
        PetscPrintf(PETSC_COMM_WORLD, "Object function value: %lf \n", ObjFunVal);

        if (ResVal_y + ResVal_u + ResVal_l < tol)
        {
            break;
        }
        VecDestroy(&dy);
        VecDestroy(&dl);
        VecDestroy(&du);
        VecDestroy(&y_new);
        VecDestroy(&y_tmp);
        VecDestroy(&u_new);
        VecDestroy(&u_tmp);
        KSPDestroy(&ksp);

        // getchar();
    }

    CleanUpPetscVar();

    MPI_Barrier(comm);
}

void DiffConv2D::CollectAndUpdateResult(Vec dx, Vec x)
{
    Vec dx_seq;
    VecScatter scatter_ctx;
    PetscReal *dx_array;
    VecScatterCreateToAll(dx, &scatter_ctx, &dx_seq);
    VecScatterBegin(scatter_ctx, dx, dx_seq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter_ctx, dx, dx_seq, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(dx_seq, &dx_array);
    MPI_Barrier(comm);

    // Update variable for visualization and residual vector computation
    //! For nonlinear iteration
    // for (uint j = 0; j < nTstep; j++)
    // {
    // 	for (uint i = 0; i < nPoint; i++)
    // 	{
    // 		int A, B;
    // 		for(uint k =0; k< state_num;k++)
    // 		{
    // 			A = i + j * nPoint;
    // 			B = i + k * nPoint + j * nPoint * state_num;
    // 			State[k][A] += PetscRealPart(dx_array[B]);
    // 		}
    // 		for(uint k =0; k< ctrl_num;k++)
    // 		{
    // 			A = i + j * nPoint;
    // 			B = i + k * nPoint + j * nPoint * ctrl_num + state_num * nPoint
    // * nTstep; 			Ctrl[k][A] += PetscRealPart(dx_array[B]);
    // 		}
    // 		for(uint k =0; k< state_num;k++)
    // 		{
    // 			A = i + j * nPoint;
    // 			B = i + k * nPoint + j * nPoint * state_num + (state_num
    // + ctrl_num) * nPoint * nTstep; 			Lambda[k][A] +=
    // PetscRealPart(dx_array[B]);
    // 		}
    // 	}
    // }
    // ! for linear system
    for (uint j = 0; j < nTstep; j++)
    {
        for (uint i = 0; i < nPoint; i++)
        {
            int A, B;
            for (uint k = 0; k < state_num; k++)
            {
                A = i + j * nPoint;
                B = i + k * nPoint + j * nPoint * state_num;
                State[k][A] = PetscRealPart(dx_array[B]);
            }
            for (uint k = 0; k < ctrl_num; k++)
            {
                A = i + j * nPoint;
                B = i + k * nPoint + j * nPoint * ctrl_num +
                    state_num * nPoint * nTstep;
                Ctrl[k][A] = PetscRealPart(dx_array[B]);
            }
            for (uint k = 0; k < state_num; k++)
            {
                A = i + j * nPoint;
                B = i + k * nPoint + j * nPoint * state_num +
                    (state_num + ctrl_num) * nPoint * nTstep;
                Lambda[k][A] = PetscRealPart(dx_array[B]);
            }
        }
    }

    VecRestoreArray(dx_seq, &dx_array);
    VecScatterDestroy(&scatter_ctx);
    VecDestroy(&dx_seq);
}

void DiffConv2D::CleanUpPetscVar()
{
    if (Y_d != NULL)
    {
        VecDestroy(&Y_d);
    }
    if (Y_bc != NULL)
    {
        VecDestroy(&Y_bc);
    }
    if (Y_ini != NULL)
    {
        VecDestroy(&Y_ini);
    }
    if (U_ini != NULL)
    {
        VecDestroy(&U_ini);
    }
    if (L_ini != NULL)
    {
        VecDestroy(&L_ini);
    }
    if (Y_k != NULL)
    {
        VecDestroy(&Y_k);
    }
    if (U_k != NULL)
    {
        VecDestroy(&U_k);
    }
    if (L_k != NULL)
    {
        VecDestroy(&L_k);
    }
    if (dX != NULL)
    {
        VecDestroy(&dX);
    }
    if (X != NULL)
    {
        VecDestroy(&X);
    }
    if (Res_nl != NULL)
    {
        VecDestroy(&Res_nl);
    }
    if (Res_Desire != NULL)
    {
        VecDestroy(&Res_Desire);
    }
    if (ResVec != NULL)
    {
        VecDestroy(&ResVec);
    }
    VecDestroy(&bsub[0]);
    VecDestroy(&bsub[1]);
    VecDestroy(&bsub[2]);

    MatDestroy(&M);
    MatDestroy(&K);
    MatDestroy(&Stable);
    MatDestroy(&Conv);
    MatDestroy(&Asubmat[0]);
    MatDestroy(&Asubmat[1]);
    MatDestroy(&Asubmat[2]);
    MatDestroy(&Asubmat[3]);
    MatDestroy(&Asubmat[4]);
    MatDestroy(&Asubmat[5]);
    MatDestroy(&Asubmat[6]);
    MatDestroy(&Asubmat[7]);
    MatDestroy(&Asubmat[8]);
    // MatDestroy(PCsubmat);
    MatDestroy(&P[0]);
    MatDestroy(&P[1]);

    // MatDestroy(&PCMat_tmp);
    // MatDestroy(&PCMat);
    MatDestroy(&TanMat_tmp);
    MatDestroy(&TanMat);
}

void DiffConv2D::VisualizeVTK_ControlMesh_Heat(const vector<Vertex2D> &spt,
                                               const vector<Element2D> &mesh,
                                               int time, int step, string fn)
{
    stringstream ss;
    ss << step << "_" << time;

    string fname = fn + "Result_CM_" + ss.str() + ".vtk";
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
                 << " " << State[0][i + time * nPoint] << "\n";
        fout << "VECTORS u float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << 0 << " " << std::setprecision(9) << 0
                 << " " << Ctrl[0][i + time * nPoint] << "\n";
        fout << "VECTORS lamda float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << 0 << " " << std::setprecision(9) << 0
                 << " " << Lambda[0][i + time * nPoint] << "\n";
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }

    // fname = fn + "Result_CM_y_"+ss.str()+".vtk";
    // fout.open(fname.c_str());
    // if (fout.is_open())
    // {
    // 	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET
    // UNSTRUCTURED_GRID\n"; 	fout << "POINTS " << spt.size() << " float\n";
    // for (i = 0; i<spt.size(); i++)
    // 	{
    // 		fout << std::setprecision(9) << spt[i].coor[0]<< std::fixed << "
    // "<< std::setprecision(9) <<spt[i].coor[1]<< std::fixed << " "<<
    // std::setprecision(9) << State[0][i + time * nPoint]<< std::fixed << "\n";
    // 	}
    // 	fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " "
    // << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
    // 	}
    // 	fout << "\nCELL_TYPES " << mesh.size() << "\n";
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "9\n";
    // 	}
    // 	fout << "POINT_DATA " << spt.size() << "\n";
    // 	fout << "VECTORS y float\n";
    // 	for (i = 0; i < spt.size(); i++)
    // 		fout << std::setprecision(9)<< State[0][i + time * nPoint] << "
    // "<< std::setprecision(9) << 0 << " " << 0 << "\n"; 	fout << "VECTORS
    // u float\n"; 	for (i = 0; i < spt.size(); i++) 		fout<< std::setprecision(9) <<
    // Ctrl[0][i+ time * nPoint] << " " << std::setprecision(9)<< 0 << " " << 0 <<
    // "\n"; 	fout << "VECTORS lamda
    // float\n"; 	for (i = 0; i < spt.size(); i++) 		fout <<
    // std::setprecision(9)<< Lambda[0][i+ time * nPoint] << " " <<
    // std::setprecision(9)<< 0 << " " << 0
    // << "\n"; 	fout.close();
    // }
    // else
    // {
    // 	cout << "Cannot open " << fname << "!\n";
    // }

    // fname = fn + "Result_CM_u_"+ss.str()+".vtk";
    // fout.open(fname.c_str());
    // if (fout.is_open())
    // {
    // 	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET
    // UNSTRUCTURED_GRID\n"; 	fout << "POINTS " << spt.size() << " float\n";
    // for (i = 0; i<spt.size(); i++)
    // 	{
    // 		fout << std::setprecision(9) << spt[i].coor[0]<< std::fixed << "
    // "<< std::setprecision(9) <<spt[i].coor[1]<< std::fixed << " "<<
    // std::setprecision(9) << Ctrl[0][i+ time * nPoint]<< std::fixed << "\n";
    // 	}
    // 	fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " "
    // << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
    // 	}
    // 	fout << "\nCELL_TYPES " << mesh.size() << "\n";
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "9\n";
    // 	}
    // 	fout << "POINT_DATA " << spt.size() << "\n";
    // 	fout << "VECTORS y float\n";
    // 	for (i = 0; i < spt.size(); i++)
    // 		fout << std::setprecision(9)<< State[0][i + time * nPoint] << "
    // "<< std::setprecision(9) << 0 << " " << 0 << "\n"; 	fout << "VECTORS
    // u float\n"; 	for (i = 0; i < spt.size(); i++) 		fout<< std::setprecision(9) <<
    // Ctrl[0][i+ time * nPoint] << " " << std::setprecision(9)<< 0 << " " << 0 <<
    // "\n"; 	fout << "VECTORS lamda
    // float\n"; 	for (i = 0; i < spt.size(); i++) 		fout <<
    // std::setprecision(9)<< Lambda[0][i+ time * nPoint] << " " <<
    // std::setprecision(9)<< 0 << " " << 0
    // << "\n"; 	fout.close();
    // }
    // else
    // {
    // 	cout << "Cannot open " << fname << "!\n";
    // }

    // fname = fn + "Result_CM_l_"+ss.str()+".vtk";
    // fout.open(fname.c_str());
    // if (fout.is_open())
    // {
    // 	fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET
    // UNSTRUCTURED_GRID\n"; 	fout << "POINTS " << spt.size() << " float\n";
    // for (i = 0; i<spt.size(); i++)
    // 	{
    // 		fout << std::setprecision(9) << spt[i].coor[0]<< std::fixed << "
    // "<< std::setprecision(9) <<spt[i].coor[1]<< std::fixed << " "<<
    // std::setprecision(9) << Lambda[0][i+ time * nPoint]<< std::fixed << "\n";
    // 	}
    // 	fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " "
    // << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
    // 	}
    // 	fout << "\nCELL_TYPES " << mesh.size() << "\n";
    // 	for (i = 0; i<mesh.size(); i++)
    // 	{
    // 		fout << "9\n";
    // 	}
    // 	fout << "POINT_DATA " << spt.size() << "\n";
    // 	fout << "VECTORS y float\n";
    // 	for (i = 0; i < spt.size(); i++)
    // 		fout << std::setprecision(9)<< State[0][i + time * nPoint] << "
    // "<< std::setprecision(9) << 0 << " " << 0 << "\n"; 	fout << "VECTORS
    // u float\n"; 	for (i = 0; i < spt.size(); i++) 		fout<< std::setprecision(9) <<
    // Ctrl[0][i+ time * nPoint] << " " << std::setprecision(9)<< 0 << " " << 0 <<
    // "\n"; 	fout << "VECTORS lamda
    // float\n"; 	for (i = 0; i < spt.size(); i++) 		fout <<
    // std::setprecision(9)<< Lambda[0][i+ time * nPoint] << " " <<
    // std::setprecision(9)<< 0 << " " << 0
    // << "\n"; 	fout.close();
    // }
    // else
    // {
    // 	cout << "Cannot open " << fname << "!\n";
    // }

    // string fname1 = fn + "velocityfield.txt";
    // ofstream fout1;
    // fout1.open(fname1.c_str());
    // if (fout1.is_open())
    // {
    // 	for (uint i = 0; i<Vel.size() / 3; i++)
    // 	{
    // 		fout1 << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 +
    // 2]
    // << "\n";
    // 	}
    // 	fout1.close();
    // }
    // else
    // {
    // 	cout << "Cannot open " << fname1 << "!\n";
    // }
}

void DiffConv2D::VisualizeVTK_ControlMesh_Burger(const vector<Vertex2D> &spt,
                                                 const vector<Element2D> &mesh,
                                                 int time, int step,
                                                 string fn)
{
    stringstream ss;
    ss << step << "_" << time;

    string fname = fn + "Result_CM_" + ss.str() + ".vtk";
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
            fout << std::setprecision(9) << State[0][i + time * nPoint] << " "
                 << std::setprecision(9) << 0 << " " << 0 << "\n";
        fout << "VECTORS u float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << Ctrl[0][i + time * nPoint] << " "
                 << std::setprecision(9) << 0 << " " << 0 << "\n";
        fout << "VECTORS lamda float\n";
        for (i = 0; i < spt.size(); i++)
            fout << std::setprecision(9) << Lambda[0][i + time * nPoint] << " "
                 << std::setprecision(9) << 0 << " " << 0 << "\n";
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }

    // string fname1 = fn + "velocityfield.txt";
    // ofstream fout1;
    // fout1.open(fname1.c_str());
    // if (fout1.is_open())
    // {
    // 	for (uint i = 0; i<Vel.size() / 3; i++)
    // 	{
    // 		fout1 << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 +
    // 2]
    // << "\n";
    // 	}
    // 	fout1.close();
    // }
    // else
    // {
    // 	cout << "Cannot open " << fname1 << "!\n";
    // }
}

void DiffConv2D::VisualizeVTK_ControlMesh(const vector<Vertex2D> &spt,
                                          const vector<Element2D> &mesh,
                                          int step, string fn)
{
    stringstream ss;
    ss << step;

    // string fname = fn + "controlmesh_VelocityPressure_"+ss.str()+".vtk";
    string fname = fn + "controlmesh_VelocityPressure.vtk";
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
            fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << spt[i].coor[2]
                 << "\n";
        }
        fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
        for (i = 0; i < mesh.size(); i++)
        {
            fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " "
                 << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << " " << mesh[i].IEN[4]
                 << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " "
                 << mesh[i].IEN[7] << '\n';
        }
        fout << "\nCELL_TYPES " << mesh.size() << '\n';
        for (i = 0; i < mesh.size(); i++)
        {
            fout << "12\n";
        }
        fout << "POINT_DATA " << Vel.size() / 3
             << "\nVECTORS VelocityField float\n";
        for (uint i = 0; i < Vel.size() / 3; i++)
        {
            fout << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 + 2]
                 << "\n";
        }
        fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
        for (uint i = 0; i < Pre.size(); i++)
        {
            fout << sqrt(Vel[i * 3] * Vel[i * 3] + Vel[i * 3 + 1] * Vel[i * 3 + 1] +
                         Vel[i * 3 + 2] * Vel[i * 3 + 2])
                 << "\n";
        }
        fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
        for (uint i = 0; i < Pre.size(); i++)
        {
            fout << Pre[i] << "\n";
        }
        fout.close();
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }

    string fname1 = fn + "velocityfield.txt";
    ofstream fout1;
    fout1.open(fname1.c_str());
    if (fout1.is_open())
    {
        for (uint i = 0; i < Vel.size() / 3; i++)
        {
            fout1 << Vel[i * 3] << " " << Vel[i * 3 + 1] << " " << Vel[i * 3 + 2]
                  << "\n";
        }
        fout1.close();
    }
    else
    {
        cout << "Cannot open " << fname1 << "!\n";
    }
}

void DiffConv2D::VisualizeVTK_PhysicalDomain(int time, int step, string fn)
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

void DiffConv2D::WriteVTK(const vector<array<double, 3>> spt,
                          const vector<double> sdisp,
                          const vector<array<int, 4>> sele, int time, int step,
                          string fn)
{
    stringstream ss;
    ss << step << "_" << time;

    string fname = fn + +"Result_physical_" + ss.str() + ".vtk";
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

void DiffConv2D::ResultCal_Bezier(double u, double v, int time,
                                  const Element2D &bzel, double pt[2],
                                  double result[result_num], double dudx[3],
                                  double &detJ)
{
    double dUdx[dim][dim];
    vector<double> Nx(bzel.IEN.size());
    vector<array<double, dim>> dNdx(bzel.IEN.size());
    vector<array<array<double, dim>, dim>> dN2dx2;
    bzel.Para2Phys(u, v, pt);
    BasisFunction(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx, dNdx, dN2dx2,
                  dUdx, detJ);
    for (uint i = 0; i < result_num; i++)
        result[i] = 0.0;

    int count = 0;
    for (uint i = 0; i < bzel.IEN.size(); i++)
    {
        for (uint j = 0; j < state_num; j++)
        {
            result[j] += Nx[i] * State[j][bzel.IEN[i] + time * nPoint];
        }
        for (uint j = 0; j < ctrl_num; j++)
        {
            result[j + state_num] += Nx[i] * Ctrl[j][bzel.IEN[i] + time * nPoint];
        }
        for (uint j = 0; j < state_num; j++)
        {
            result[j + state_num + ctrl_num] +=
                Nx[i] * Lambda[j][bzel.IEN[i] + time * nPoint];
        }
    }
}

void DiffConv2D::DebugSubMat()
{
    IS isg[2];

    ISCreateStride(PETSC_COMM_WORLD, (state_num + ctrl_num) * nPoint * nTstep, 0,
                   1, &isg[0]);
    ISCreateStride(PETSC_COMM_WORLD, state_num * nPoint * nTstep,
                   (state_num + ctrl_num) * nPoint * nTstep, 1, &isg[1]);
    Mat Asub[4];
    MatCreateSubMatrix(TanMat, isg[0], isg[0], MAT_INITIAL_MATRIX, &Asub[0]);
    MatCreateSubMatrix(TanMat, isg[0], isg[1], MAT_INITIAL_MATRIX, &Asub[1]);
    MatCreateSubMatrix(TanMat, isg[1], isg[0], MAT_INITIAL_MATRIX, &Asub[2]);
    MatCreateSubMatrix(TanMat, isg[1], isg[1], MAT_INITIAL_MATRIX, &Asub[3]);

    ierr = PetscObjectSetName((PetscObject)Asub[0], "TanMat00");
    ierr = PetscObjectSetName((PetscObject)Asub[1], "TanMat01");
    ierr = PetscObjectSetName((PetscObject)Asub[2], "TanMat10");
    ierr = PetscObjectSetName((PetscObject)Asub[3], "TanMat11");

    DebugVisualizeMat(Asub[0], work_dir + "debug/TanMat00_afterBC.m");
    DebugVisualizeMat(Asub[1], work_dir + "debug/TanMat01_afterBC.m");
    DebugVisualizeMat(Asub[2], work_dir + "debug/TanMat10_afterBC.m");
    DebugVisualizeMat(Asub[3], work_dir + "debug/TanMat11_afterBC.m");
}

// For 3D
// void DiffConv2D::BasisFunction(double u, double v, double w, int nen, const
// vector<array<double, 3>>& pt, const vector<array<double, 64>> &cmat,
// vector<double> &Nx, vector<array<double, 3>> &dNdx,
// vector<array<array<double, 3>, 3>> &dN2dx2, double dudx[3][3], double& detJ)
// {
// 	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. -
// u)*u, 3.*(1. - u)*u*u, u*u*u }; 	double Nv[4] = { (1. - v)*(1. - v)*(1. -
// v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v }; 	double Nw[4] = {
// (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };

// 	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2.
// - 3.*u)*u, 3.*u*u }; 	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3.
// - 12.*v
// + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v }; 	double dNdw[4] = { -3.*(1. -
// w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };

// 	double dN2du2[4] = { 6.*(1. - u),-12. + 18.*u,6. - 18.*u,6.*u };
// 	double dN2dv2[4] = { 6.*(1. - v),-12. + 18.*v,6. - 18.*v,6.*v };
// 	double dN2dw2[4] = { 6.*(1. - w),-12. + 18.*w,6. - 18.*w,6.*w };

// 	double dNdt[bzpt_num][3];
// 	double dN2dt2[bzpt_num][3][3];
// 	double Nx_bz[bzpt_num];
// 	double dNdx_bz[bzpt_num][3];
// 	double dN2dx2_bz[bzpt_num][3][3];

// 	Nx.clear();
// 	dNdx.clear();
// 	dN2dx2.clear();
// 	Nx.resize(nen, 0);
// 	dNdx.resize(nen, { 0 });
// 	dN2dx2.resize(nen, { {0} });

// 	int i, j, k, a, b, c, loc;
// 	loc = 0;
// 	for (i = 0; i<4; i++){
// 		for (j = 0; j<4; j++){
// 			for (k = 0; k < 4; k++)	{
// 				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
// 				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
// 				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
// 				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
// 				dN2dt2[loc][0][0] = dN2du2[k] * Nv[j] * Nw[i];
// dN2dt2[loc][0][1] = dNdu[k] * dNdv[j] * Nw[i];
// dN2dt2[loc][0][2] = dNdu[k] * Nv[j] * dNdw[i];
// dN2dt2[loc][1][0] = dNdu[k] * dNdv[j] * Nw[i];
// dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j] * Nw[i];
// dN2dt2[loc][1][2] = Nu[k] * dNdv[j] * dNdw[i];
// dN2dt2[loc][2][0] = dNdu[k] * Nv[j] * dNdw[i]; dN2dt2[loc][2][1] = Nu[k] *
// dNdv[j] * dNdw[i]; dN2dt2[loc][2][2] = Nu[k] * Nv[j] * dN2dw2[i];
// loc++;
// 			}
// 		}
// 	}

// 	double dxdt[3][3] = { {0} };
// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (a = 0; a<3; a++)
// 			for (b = 0; b<3; b++)
//  				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];

// 	double dtdx[3][3] = { { 0 } };
// 	Matrix3DInverse(dxdt, dtdx);

// 	//1st derivatives
// 	for (i = 0; i<bzpt_num; i++){
// 		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0]
// +
// dNdt[i][2] * dtdx[2][0]; 		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] +
// dNdt[i][1]
// * dtdx[1][1] + dNdt[i][2] * dtdx[2][1]; 		dNdx_bz[i][2] = dNdt[i][0]
// * dtdx[0][2] + dNdt[i][1] * dtdx[1][2] + dNdt[i][2] * dtdx[2][2];
// 	}
// 	for (int i = 0; i < 3; i++)
// 		for (int j = 0; j < 3; j++)
// 			dudx[i][j] = dtdx[i][j];

// 	detJ = MatrixDet(dxdt);
// 	detJ = 0.125*detJ;

// 	//2nd derivatives
// 	double dx2dt2[3][9] = { {0} };
// 	double dt2dx2[3][9] = { {0} };
// 	for (int l = 0; l < 3; l++)
// 		for (loc = 0; loc<bzpt_num; loc++)
// 			for (a = 0; a<3; a++)
// 				for (b = 0; b<3; b++)
// 					dx2dt2[l][3*a+b] += pt[loc][l] *
// dN2dt2[loc][a][b];

// 	for (i = 0; i < 3; i++)
// 		for (j = 0; j < 3; j++)
// 			for (k = 0; k < 3; k++)
// 				for (a = 0; a < 3; a++)
// 					for (b = 0; b < 3; b++)
// 						for (c = 0; c < 3; c++)
// 							dt2dx2[c][3*i+j] -=
// dx2dt2[k][3*a+b]*dtdx[a][i]*dtdx[b][j]*dtdx[c][k];

// 	for (loc = 0; loc < bzpt_num; loc++)
// 		for (i = 0; i < 3; i++)
// 			for (j = 0; j < 3; j++)
// 				dN2dx2_bz[loc][i][j] = 0.;

// 	for (loc = 0; loc < bzpt_num; loc++){
// 		for (i = 0; i < 3; i++)	{
// 			for (j = 0; j < 3; j++)	{
// 				for (a = 0; a<3; a++){
// 					for (b = 0; b<3; b++){
// 						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b]
// * dtdx[a][i]*dtdx[b][j];
// 					}
// 					dN2dx2_bz[loc][i][j] += dNdt[loc][a] *
// dt2dx2[a][3*i+j];
// 				}
// 			}
// 		}
// 	}

// 	for (i = 0; i < nen; i++){
// 		for (j = 0; j < bzpt_num; j++)	{
// 			Nx[i] += cmat[i][j] * Nx_bz[j];
// 			for (int m = 0; m < 3; m++)	{
// 				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
// 				for (int n = 0; n < 3; n++)	{
// 					dN2dx2[i][m][n] += cmat[i][j] *
// dN2dx2_bz[j][m][n];
// 				}
// 			}
// 		}
// 	}
// }

// void DiffConv2D::VisualizeVTK_PhysicalDomain(int step, string fn)
// {
// 	vector<array<double, 3>> spt_all;//sample points
// 	vector<double> sresult_all;
// 	vector<array<int, 8>> sele_all;
// 	double detJ;
// 	int num_bzmesh_ele = bzmesh_process.size();
// 	double spt_proc[num_bzmesh_ele * 24];
// 	double sresult_proc[num_bzmesh_ele * 32];
// 	int sele_proc[num_bzmesh_ele * 8];
// 	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
// 	{
// 		int ns(2);
// 		vector<double> su(ns);
// 		for (int i = 0; i < ns; i++)
// 		{
// 			su[i] = double(i) / (double(ns) - 1.);
// 		}

// 		int loc(0);
// 		for (int a = 0; a<ns; a++)
// 		{
// 			for (int b = 0; b<ns; b++)
// 			{

// 					double pt1[3], dudx[3];
// 					double result[4];
// 					ResultCal_Bezier(su[b], su[a],
// bzmesh_process[e],
// pt1, result, dudx, detJ); 					spt_proc[24 * e + loc*3 +
// 0] = pt1[0]; 					spt_proc[24
// * e + loc*3 + 1] = pt1[1]; 					spt_proc[24 * e + loc*3 +
// 2] = pt1[2]; 					sresult_proc[32 * e + loc * 4 + 0] = result[0];
// sresult_proc[32 * e + loc * 4 + 1] = result[1];
// sresult_proc[32 * e + loc * 4 + 2] = result[2];
// sresult_proc[32 * e + loc * 4 + 3] = result[3]; 					loc++;

// 			}
// 		}
// 		int nns[2] = { ns*ns*ns, ns*ns };
// 		for (int a = 0; a<ns - 1; a++)
// 		{
// 			for (int b = 0; b<ns - 1; b++)
// 			{
// 				for (int c = 0; c < ns - 1; c++)
// 				{
// 					sele_proc[8 * e + 0] = 8 * e + a*nns[1] + b*ns
// +
// c; 					sele_proc[8 * e + 1] = 8 * e + a*nns[1] + b*ns +
// c + 1; 					sele_proc[8 * e + 2] = 8 * e +
// a*nns[1] + (b + 1)*ns + c + 1; 					sele_proc[8 * e + 3]
// = 8 * e +
// a*nns[1] + (b + 1)*ns + c; 					sele_proc[8 * e + 4] = 8
// * e + (a + 1)*nns[1] +
// b*ns + c; 					sele_proc[8 * e + 5] = 8 * e + (a +
// 1)*nns[1] + b*ns + c + 1; 					sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b
// +
// 1)*ns + c + 1; 					sele_proc[8 * e + 7] = 8 * e + (a +
// 1)*nns[1] + (b + 1)*ns + c;
// 				}
// 			}
// 		}
// 	}

// 	double *spts= NULL;
// 	double *sresults=NULL;
// 	int *seles=NULL;
// 	int *displs_spts = NULL;
// 	int *displs_sresults = NULL;
// 	int *displs_seles = NULL;
// 	int *num_bzmesh_eles=NULL;
// 	int *recvcounts_spts = NULL;
// 	int *recvcounts_sresults = NULL;
// 	int *recvcounts_seles = NULL;

// 	if (comRank == 0)
// 	{
// 		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
// 		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
// 		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
// 		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
// 	}
// 	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0,
// PETSC_COMM_WORLD); 	MPI_Barrier(comm);

// 	if (comRank == 0)
// 	{
// 		spts = (double*)malloc(sizeof(double) * 24 * n_bzmesh);
// 		sresults = (double*)malloc(sizeof(double) * 32 * n_bzmesh);
// 		seles = (int*)malloc(sizeof(int) * 8 * n_bzmesh);

// 		displs_spts = (int*)malloc(nProcess * sizeof(int));
// 		displs_sresults = (int*)malloc(nProcess * sizeof(int));
// 		displs_seles = (int*)malloc(nProcess * sizeof(int));
// 		displs_spts[0] = 0;
// 		displs_sresults[0] = 0;
// 		displs_seles[0] = 0;

// 		for (int i = 1; i<nProcess; i++) {
// 			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1]
// * 24; 			displs_sresults[i] = displs_sresults[i - 1] +
// num_bzmesh_eles[i - 1] *
// 32; 			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i
// - 1] * 8;
// 		}

// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			recvcounts_spts[i] = num_bzmesh_eles[i] * 24;
// 			recvcounts_sresults[i] = num_bzmesh_eles[i] * 32;
// 			recvcounts_seles[i] = num_bzmesh_eles[i] * 8;
// 		}
// 	}

// 	MPI_Gatherv(spt_proc, num_bzmesh_ele * 8 * 3, MPI_DOUBLE, spts,
// recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 8 * 4, MPI_DOUBLE, sresults,
// recvcounts_sresults, displs_sresults ,MPI_DOUBLE, 0, PETSC_COMM_WORLD);
// 	MPI_Gatherv(sele_proc, num_bzmesh_ele * 8, MPI_INT, seles,
// recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);

// 	if (comRank == 0)
// 	{
// 		for (int i = 0; i < n_bzmesh; i++)
// 		{
// 			for (int j = 0; j < 8; j++)
// 			{
// 				array<double, 3> pt = { spts[i * 24 + j * 3 + 0],
// spts[i
// * 24 + j * 3 + 1], spts[i * 24 + j * 3 + 2] };
// spt_all.push_back(pt); 				sresult_all.push_back(sresults[i * 32 + j * 4 + 0]);
// sresult_all.push_back(sresults[i * 32 + j * 4 + 1]);
// 				sresult_all.push_back(sresults[i * 32 + j * 4 +
// 2]); 				sresult_all.push_back(sresults[i * 32 + j * 4 +
// 3]);
// 			}

// 		}
// 		int sum_ele = 0;
// 		int pstart = 0;
// 		for (int i = 0; i < nProcess; i++)
// 		{
// 			for (int e = 0; e < num_bzmesh_eles[i]; e++)
// 			{
// 				array<int, 8> el;
// 				el[0] = pstart + seles[8 * sum_ele + 0];
// 				el[1] = pstart + seles[8 * sum_ele + 1];
// 				el[2] = pstart + seles[8 * sum_ele + 2];
// 				el[3] = pstart + seles[8 * sum_ele + 3];
// 				el[4] = pstart + seles[8 * sum_ele + 4];
// 				el[5] = pstart + seles[8 * sum_ele + 5];
// 				el[6] = pstart + seles[8 * sum_ele + 6];
// 				el[7] = pstart + seles[8 * sum_ele + 7];
// 				sele_all.push_back(el);
// 				sum_ele++;
// 			}
// 			pstart = pstart + num_bzmesh_eles[i] * 8;
// 		}
// 		cout << "Visualizing in Physical Domain...\n";
// 		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
// 	}

// }

// void DiffConv2D::WriteVTK(const vector<array<double, 3>> spt, const
// vector<double> sdisp, const vector<array<int,8>> sele, int step, string fn)
// {
// 	string fname = fn + "_VelocityPressure.vtk";
// 	ofstream fout;
// 	fout.open(fname.c_str());
// 	unsigned int i;
// 	if (fout.is_open())
// 	{
// 		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET
// UNSTRUCTURED_GRID\n"; 		fout << "POINTS " << spt.size() << "
// float\n"; 		for (i = 0; i<spt.size(); i++)
// 		{
// 			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2]
// <<
// "\n";
// 		}
// 		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() <<
// '\n'; 		for (i = 0; i<sele.size(); i++)
// 		{
// 			fout << "8 " << sele[i][0] << " " << sele[i][1] << " "
// << sele[i][2] << " " << sele[i][3]
// 				<< " " << sele[i][4] << " " << sele[i][5] << " "
// << sele[i][6] << " " << sele[i][7] << '\n';
// 		}
// 		fout << "\nCELL_TYPES " << sele.size() << '\n';
// 		for (i = 0; i<sele.size(); i++)
// 		{
// 			fout << "12\n";
// 		}
// 		fout << "POINT_DATA " << sdisp.size() / 4 << "\nVECTORS
// VelocityField float\n"; 		for (uint i = 0; i<sdisp.size() / 4;
// i++)
// 		{
// 			fout << sdisp[i * 4] << " " << sdisp[i * 4 + 1] << " "
// << sdisp[i * 4 + 2] << "\n";
// 		}
// 		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE
// default\n"; 		for (uint i = 0; i<sdisp.size() / 4; i++)
// 		{
// 			fout << sqrt(sdisp[i * 4] * sdisp[i * 4] + sdisp[i * 4 + 1]
// * sdisp[i * 4 + 1] + sdisp[i * 4 + 2] * sdisp[i * 4 + 2]) << "\n";
// 		}
// 		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
// 		for (uint i = 0; i<sdisp.size() / 4; i++)
// 		{
// 			fout << sdisp[4 * i + 3] << "\n";
// 		}
// 		fout.close();
// 	}
// 	else
// 	{
// 		cout << "Cannot open " << fname << "!\n";
// 	}
// }

// void DiffConv2D::ResultCal_Bezier(double u, double v, const Element2D& bzel,
// double pt[3], double result[4], double dudx[3], double& detJ)
// {
// 	double dUdx[dim][dim];
// 	vector<double> Nx(bzel.IEN.size());
// 	vector<array<double, dim>> dNdx(bzel.IEN.size());
// 	vector<array<array<double, dim>, dim>> dN2dx2;
// 	bzel.Para2Phys(u, v, pt);
// 	BasisFunction(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx,
// dNdx,dN2dx2, dUdx, detJ); 	result[0] = 0.; result[1] = 0.; result[2] = 0.;
// result[3] = 0.; 	for (uint i = 0; i < bzel.IEN.size(); i++)	{
// result[0] += Nx[i] * (Vel[dim* bzel.IEN[i] + 0]); 		result[1] +=
// Nx[i] * (Vel[dim*
// bzel.IEN[i] + 1]); 		result[2] += Nx[i] * (Vel[dim* bzel.IEN[i] +
// 2]); 		result[3] += Nx[i] * (Pre[bzel.IEN[i]]);
// 	}
// }
// }
