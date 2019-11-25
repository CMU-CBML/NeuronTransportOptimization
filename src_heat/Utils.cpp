#include "Utils.h"

double MatrixDet(double dxdt[2][2])
{
    double det = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
    return det;
}

double MatrixDet(double dxdt[3][3])
{
    double det = dxdt[0][0] * dxdt[1][1] * dxdt[2][2] + dxdt[0][1] * dxdt[1][2] * dxdt[2][0] + dxdt[0][2] * dxdt[2][1] * dxdt[1][0] -
                 (dxdt[0][2] * dxdt[1][1] * dxdt[2][0] + dxdt[0][0] * dxdt[1][2] * dxdt[2][1] + dxdt[1][0] * dxdt[0][1] * dxdt[2][2]);

    return det;
}

void Matrix2DInverse(double dxdt[2][2], double dtdx[2][2])
{
    double det = MatrixDet(dxdt);
    dtdx[0][0] = 1.0 / det * (dxdt[1][1]);
    dtdx[0][1] = 1.0 / det * (-dxdt[0][1]);
    dtdx[1][0] = 1.0 / det * (-dxdt[1][0]);
    dtdx[1][1] = 1.0 / det * (dxdt[0][0]);
}

void Matrix3DInverse(double dxdt[3][3], double dtdx[3][3])
{
    double det = MatrixDet(dxdt);
    dtdx[0][0] = 1 / det * (dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]);
    dtdx[0][1] = 1 / det * (dxdt[2][1] * dxdt[0][2] - dxdt[0][1] * dxdt[2][2]);
    dtdx[0][2] = 1 / det * (dxdt[0][1] * dxdt[1][2] - dxdt[1][1] * dxdt[0][2]);
    dtdx[1][0] = 1 / det * (dxdt[2][0] * dxdt[1][2] - dxdt[1][0] * dxdt[2][2]);
    dtdx[1][1] = 1 / det * (dxdt[0][0] * dxdt[2][2] - dxdt[0][2] * dxdt[2][0]);
    dtdx[1][2] = 1 / det * (dxdt[1][0] * dxdt[0][2] - dxdt[0][0] * dxdt[1][2]);
    dtdx[2][0] = 1 / det * (dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
    dtdx[2][1] = 1 / det * (dxdt[0][1] * dxdt[2][0] - dxdt[0][0] * dxdt[2][1]);
    dtdx[2][2] = 1 / det * (dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0]);
}

void DebugMatVector(vector<vector<double>> mat_debug, int e, string fname)
{
    int i, j;
    ofstream fout;

    fout.open(fname.c_str(), ios::app);

    if (fout.is_open())
    {
        fout << "Element Index: " << e << "\n";
        fout << "Size: " << mat_debug.size() << " " << mat_debug.size() << "\n";
        for (i = 0; i < mat_debug.size(); i++)
        {
            for (j = 0; j < mat_debug[i].size(); j++)
            {
                fout << mat_debug[i][j] << " ";
            }
            fout << "\n";
        }
    }
    else
    {
        cout << "Cannot open " << fname << "!\n";
    }
}

void DebugVisualizeMat(Mat mat_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = MatView(mat_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeVec(Vec vec_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = VecView(vec_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeIS(IS is_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = ISView(is_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeMatSelf(Mat mat_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = MatView(mat_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeVecSelf(Vec vec_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = VecView(vec_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeISSelf(IS is_debug, string fname)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, fname.c_str(), &viewer);
    ierr = PetscViewerPushFormat(viewer, format);
    ierr = ISView(is_debug, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = PetscViewerDestroy(&viewer);
}

