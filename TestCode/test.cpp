#include <petsc.h>
#include <petscksp.h>
 
#include "petscsys.h"   
#include "petscmat.h"
#include <petscis.h>
#include <petscviewer.h>
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
static char help[] = "Demonstrates creating a stride index set.\n\n";

const int Npoint = 9;
const int Ntime = 5;
const int Nctrl = 4;
const int Nstat = 7;
const int Npena = 7;

PetscErrorCode GetPosition(PetscInt row, PetscInt n_point, PetscInt n_var, PetscInt n_tstep, PetscInt &i_point, PetscInt &i_var, PetscInt &i_tstep)
{
    PetscInt n;
    PetscFunctionBeginUser;
    /* cell number n=j*nx+i has position (i,j) in grid */
    i_tstep = row / (n_var * n_point);
    i_var = row % (n_var * n_point) / n_point;
    i_point = row % (n_var * n_point) % n_point;
    PetscFunctionReturn(0);
}

static PetscErrorCode FormSubA(Mat &A, PetscInt nx, PetscInt ny,  PetscInt flag)
{
    PetscInt ind_point, ind_var, ind_time;
    PetscInt n_row = Npoint * Ntime * nx;
    PetscInt n_col = Npoint * Ntime * ny;
    PetscInt row, col, start, end;   

    MatCreate(PETSC_COMM_WORLD, &A);
    // MatSetOptionsPrefix(A,"A00_");
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, n_row, n_col);
    MatSetType(A,MATMPIAIJ);
    MatMPIAIJSetPreallocation(A, Npoint * ny, NULL, Npoint * ny ,NULL);
    MatGetOwnershipRange(A,&start,&end);
    
    if(flag !=0)
    {
        for(row = start; row < end; row++){
        GetPosition(row, Npoint, Nstat, Ntime, ind_point, ind_var, ind_time);
        PetscInt col_start = ind_time * Npoint * ny;
        PetscInt col_end = (ind_time+1) * Npoint * ny;
        for(col = col_start; col < col_end; col++)
        {
            PetscScalar    vals[] = {ind_point*1.0 + ind_var*10.0 + ind_time*100.0};
            MatSetValues(A, 1, &row, 1, &col, vals, INSERT_VALUES);
        }
    }
    }
    

    // for(int i=0;i<n_row;i++){
    //     for(int j=0;j<n_col;j++){
    //         row = i;
    //         col = j;
    //         PetscScalar vals[] = {i+j};
    //         MatSetValues(A,1,&row,1,&col,vals,INSERT_VALUES);
    //     }
    // }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0);
}

// static PetscErrorCode FormSubA00(Mat A, PetscInt nx, PetscInt ny)
// {
//     PetscInt ind_point, ind_var, ind_time;
//     PetscInt n_row = Npoint * Ntime * nx;
//     PetscInt n_col = Npoint * Ntime * ny;
//     PetscInt row, col, start, end;   

//     MatCreate(PETSC_COMM_WORLD, &A);
//     // MatSetOptionsPrefix(A,"A00_");
//     MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, n_row, n_col);
//     MatSetType(A,MATMPIAIJ);
//     MatMPIAIJSetPreallocation(A, Npoint * ny, NULL, Npoint * ny ,NULL);
//     MatGetOwnershipRange(A,&start,&end);
    
//     for(row = start; row < end; row++){
//         GetPosition(row, Npoint, Nstat, Ntime, ind_point, ind_var, ind_time);
//         PetscInt col_start = ind_time * Npoint * ny;
//         PetscInt col_end = (ind_time+1) * Npoint * ny;
//         for(col = col_start; col < col_end; col++)
//         {
//             PetscScalar    vals[] = {ind_point*1.0 + ind_var*10.0 + ind_time*100.0};
//             MatSetValues(A, 1, &row, 1, &col, vals, INSERT_VALUES);
//         }
//     }

//     // for(int i=0;i<n_row;i++){
//     //     for(int j=0;j<n_col;j++){
//     //         row = i;
//     //         col = j;
//     //         PetscScalar vals[] = {i+j};
//     //         MatSetValues(A,1,&row,1,&col,vals,INSERT_VALUES);
//     //     }
//     // }
//     MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//     return(0);
// }

int main(int argc, char **argv)
{
    PetscInt       i,n,first,step;
    PetscErrorCode ierr;
    IS             set1, set2, set3;
    Mat A; 
    Mat B;
    Mat a[9];
    Mat b[4];
    const PetscInt *indices;

    ierr = PetscInitialize(&argc,&argv,(char*)0,help);
    if (ierr) return ierr;
    
    //Create 3 IS sets
    n     = 7;
    first = 0;
    step  = 1;

    ISCreateStride(PETSC_COMM_WORLD,n,first,step,&set1);

    n     = 4;
    first = 7;
    step  = 1;

    ISCreateStride(PETSC_COMM_WORLD,n,first,step,&set2);

    n     = 7;
    first = 11;
    step  = 1;

    ISCreateStride(PETSC_COMM_WORLD,n,first,step,&set3);
   
    IS is_row[]={set1, set2, set3};
    IS is_col[]={set1, set2, set3};
    // ISLocalToGlobalMapping rmapping, cmapping;

    // ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, );


    // ierr = MatCreate(PETSC_COMM_WORLD, &A); 
	// ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 18, 18);
	// ierr = MatSetType(A, MATMPIAIJ);
	// ierr = MatMPIAIJSetPreallocation(A, 18, NULL, 18, NULL);

    ierr = FormSubA(B,Nstat,Nstat,1);

    PetscPrintf(PETSC_COMM_WORLD,  "finish Matrix B\n");
    PetscViewer viewer;
    PetscViewerFormat format =  PETSC_VIEWER_ASCII_MATLAB;
    

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",&viewer);
    ierr = PetscViewerPushFormat(viewer,format);
    ierr = MatView(B,viewer);

    ierr = FormSubA(a[0],Nstat,Nstat,1);
    ierr = FormSubA(a[1],Nstat,Nctrl,0);
    ierr = FormSubA(a[2],Nstat,Npena,1);
    ierr = FormSubA(a[3],Nctrl,Nstat,0);
    ierr = FormSubA(a[4],Nctrl,Nctrl,1);
    ierr = FormSubA(a[5],Nctrl,Npena,1);
    ierr = FormSubA(a[6],Npena,Nstat,1);
    ierr = FormSubA(a[7],Npena,Nctrl,1);
    ierr = FormSubA(a[8],Npena,Npena,0);
    ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, a, &A);
    PetscPrintf(PETSC_COMM_WORLD,  "finish Matrix A\n");
    // PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matA.output",&viewer);
    ierr = MatView(A,viewer);

    ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, a, &b[0]);
    ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, a, &b[1]);
    ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, a, &b[2]);
    ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, a, &b[3]);
    ierr = MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, b, &B);
    PetscPrintf(PETSC_COMM_WORLD,  "finish Matrix B\n");
    // PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matB.output",&viewer);
    ierr = MatView(B,viewer);


    // Mat C;
    // Mat c[2];
    // IS is_r, is_c;
    // // PetscInt id_r[];
    // // PetscInt id_c[];

    // ierr = MatCreateNest(PETSC_COMM_WORLD, 3, NULL, 3, NULL, NULL, &C);

    // PetscInt id_r1[] = {0};
    // PetscInt id_c1[] = {0, 2};
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 1, id_r1, PETSC_COPY_VALUES, &is_r);
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 2, id_c1, PETSC_COPY_VALUES, &is_c);
    // ierr = FormSubA(c[0],Nstat,Nstat,1);
    // ierr = FormSubA(c[1],Nstat,Npena,1);
    // MatNestSetSubMats(C,1,&is_r,2,&is_c,c);


    // PetscInt id_r2[] = {1};
    // PetscInt id_c2[] = {1, 2};
    // ierr = FormSubA(c[0],Nctrl,Nctrl,1);
    // ierr = FormSubA(c[1],Nctrl,Npena,1);
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 1, id_r2, PETSC_COPY_VALUES, &is_r);
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 2, id_c2, PETSC_COPY_VALUES, &is_c);
    // MatNestSetSubMats(C,1,&is_r,2,&is_c,c);

    // PetscInt id_r3[] = {2};
    // PetscInt id_c3[] = {0, 1};
    // ierr = FormSubA(c[0],Npena,Nstat,1);
    // ierr = FormSubA(c[1],Npena,Nctrl,1);
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 1, id_r3, PETSC_COPY_VALUES, &is_r);
    // ierr = ISCreateGeneral(PETSC_COMM_WORLD, 2, id_c3, PETSC_COPY_VALUES, &is_c);
    // MatNestSetSubMats(C,1,&is_r,2,&is_c,c);
    
    // // ierr = MatCreateNest(PETSC_COMM_WORLD, 3, &is_r, 3, &is_c, c, &C);

 

    // PetscPrintf(PETSC_COMM_WORLD,  "finish Matrix C\n");
    // // PetscViewer viewer;
    // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matC.output",&viewer);
    // ierr = MatView(C,viewer);


    // ierr = MatView(A, 	PETSC_VIEWER_STDOUT_WORLD);
    // ISGetIndices(set1,&indices);
    // PetscPrintf(PETSC_COMM_WORLD,"Printing indices directly\n");
    // for (i=0; i<n; i++) {
    //     // cout << indices[i] << endl;
    //   PetscPrintf(PETSC_COMM_WORLD,"%d\n",indices[i]);
    // }
    // ISRestoreIndices(set1,&indices);

    // ISGetIndices(set2,&indices);
    // PetscPrintf(PETSC_COMM_WORLD,"Printing indices directly\n");
    // for (i=0; i<4; i++) {
    //     // cout << indices[i] << endl;
    //   PetscPrintf(PETSC_COMM_WORLD,"%d\n",indices[i]);
    // }
    // ISRestoreIndices(set2,&indices);

    // ISGetIndices(set3,&indices);
    // PetscPrintf(PETSC_COMM_WORLD,"Printing indices directly\n");
    // for (i=0; i<n; i++) {
    //     // cout << indices[i] << endl;
    //   PetscPrintf(PETSC_COMM_WORLD,"%d\n",indices[i]);
    // }
    // ISRestoreIndices(set3,&indices);
  



    ISDestroy(&set1);
    ISDestroy(&set2);
    ISDestroy(&set3);
    PetscFinalize();
    return ierr;
}

// int main(int argc, char **argv)
// {
//     PetscInt       i,n,first,step;
//     PetscErrorCode ierr;
//     IS             set;
//     const PetscInt *indices;

//     ierr = PetscInitialize(&argc,&argv,(char*)0,help);
//     if (ierr) return ierr;
    
//     n     = 10;
//     first = 3;
//     step  = 2;

//     ISCreateStride(PETSC_COMM_WORLD,n,first,step,&set);
//     // ISView(set,PETSC_VIEWER_STDOUT_SELF);


//     ISGetIndices(set,&indices);
//     PetscPrintf(PETSC_COMM_WORLD,"Printing indices directly\n");
//     for (i=0; i<n; i++) {
//         // cout << indices[i] << endl;
//       PetscPrintf(PETSC_COMM_WORLD,"%d\n",indices[i]);
//     }
//     ISRestoreIndices(set,&indices);
//     /*
//         Determine information on stride
//     */
//     ISStrideGetInfo(set,&first,&step);
//     if (first != 3 || step != 2) SETERRQ(PETSC_COMM_SELF,1,"Stride info not correct!\n");
//     ISDestroy(&set);
//     PetscFinalize();
//     return ierr;
// }