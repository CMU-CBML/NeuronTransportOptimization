#include "UserSetting2D.h"

UserSetting2D::UserSetting2D()
{
	comm = MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;
}

UserSetting2D::~UserSetting2D()
{
}

bool UserSetting2D::GetReadDesire() const
{
	return (ReadDesire == true) ? true : false;
}

void UserSetting2D::InitializeUserSetting(string path_in, string path_out)
{
	int rank, nProcs;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

	ele_process.resize(nProcs);

	work_dir = path_in;
	output_dir = path_out;

	PetscPrintf(PETSC_COMM_WORLD, "Working Directory: %s\n", output_dir.c_str());
	PetscPrintf(PETSC_COMM_WORLD, "nProcs: %d\n", nProcs);

	string fn_mesh(path_in + "controlmesh.vtk");
	string fn_bzmeshinfo(path_in + "bzmeshinfo.txt");
	string fn_bz(path_in + "bzmeshinfo.txt.epart." + to_string(nProcs));

	string fn_IC(path_out + "InitialCondition.txt");
	string fn_IC_visual(path_out + "InitialCondition.vtk");
	string fn_BC(path_out + "BoundaryCondition.txt");
	string fn_BC_visual(path_out + "BoundaryCondition.vtk");
	string fn_Desire(path_out + "DesireState.txt");
	string fn_Desire_visual(path_out + "DesireState.vtk");

	string fn_velocity(path_in + "initial_velocityfield.txt");
	string fn_parameter(path_in + "simulation_parameter.txt");

	SetVariables(fn_parameter);

	// SetInitialCondition(int ndof, vector<double> &Vel0, vector<double> &Pre0, vector<Vertex2D> &pts, const vector<array<double, 2>> velocity_node);
	ReadMesh(fn_mesh);
	Readbzmeshinfo(fn_bzmeshinfo);
	// CreateDMPlexMesh(fn_bzmeshinfo);

	// ReadVelocityField(fn_velocity);
	AssignProcessor(fn_bz);

	// RunL2Projection();

	SetInitialCondition(fn_IC, fn_IC_visual);
	SetDesireState(fn_Desire, fn_Desire_visual);
	SetBoundaryCondition(fn_BC, fn_BC_visual);
}

void UserSetting2D::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";

	// cout << fname_cmat <<endl;

	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open())
	{
		fin_cmat >> neles;
		// cout << neles << endl;
		bzmesh_process.resize(ele_process[comRank].size());
		for (int i = 0; i < neles; i++)
		{
			if (i == ele_process[comRank][add])
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

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	// cout << fname_bzpt <<endl;
	add = 0;
	if (fin_bzpt.is_open())
	{
		fin_bzpt >> npts;
		// cout << npts << endl;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++)
		{
			if (e == ele_process[comRank][add])
			{
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++)
				{
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> bzmesh_process[add].pts[i][2];
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
}

void UserSetting2D::DesireStateFunction(double x, double y, double z, double t, double result[state_num]) const
{
	for (int i = 0; i < state_num; i++)
		result[i] = 0;

	// ! Bifurcation desire function
	double tt = sqrt(x * x + y * y);
	double tx, ty;
	double theta = atan2(y, x) * 180.0 / PI;
	double vmax;

	if (theta >= 120 || theta <= -120)
	{
		result[0] = var[3];
	}
	if (theta > 0 && theta < 120)
	{
		result[0] = var[3] * (1.0 / 6.0 * tt * tt - 2.0 / 3.0 * tt + 1);
	}
	if (theta <= 0 && theta > -120)
	{
		result[0] = var[3];
	}

	if (theta >= 120 || theta <= -120)
	{
		tx = -x;
		ty = y;

		vmax = var[3] * (tx + 4) / 6;
		// vmax = tx/3.0+2.0;
		result[0] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
	}
	if (theta > 0 && theta < 120)
	{
		tx = x * 0.5 + y * 0.5 * sqrt(3);
		ty = -x * 0.5 * sqrt(3) + y * 0.5;

		vmax = var[3] * (-tx + 5) / 10;

		// vmax = 3.0 * (-tx + 3.1)/6 + 0.2;
		// vmax = 1.0/6.0*(tx-1.0) *(tx-3.0) + 1.0;
		result[0] = 0.9 * vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		if (abs(tx - 2.0) < 0.2)
			// result[0] = 0.4 * vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
			result[0] = 1.5 * vmax;
	}
	if (theta <= 0 && theta > -120)
	{

		tx = x * 0.5 - y * 0.5 * sqrt(3);
		ty = x * 0.5 * sqrt(3) + y * 0.5;
		vmax = 3.0 * (-tx + 3.1) / 6 + 0.2;
		vmax = 1.4 * var[3] * (-tx + 5) / 10;
		result[0] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
	}

	// ! One pipe desire function
	// // result[0] = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0; // concentration
	// // 													  // result = 3.0 / 50.0 * x0 * x0 - 4.0 / 5.0 * x0 + 1.0; // velocity
	// // 													  // result[0] = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0;
	// // 													  // result[0] = exp(-x / var[0]) - 0.1 * x + 2.0;
	// double vmax = 0.0, vx = 0.0, vy = 0.0;
	// // vmax = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0;
	// // vmax = 1.0 / 70.0 * (x - 12) * (x - 12) + 33 / 35 + 0.2;
	// // // vmax = 1.0 + (10 - x) / 5;
	// vmax = 2.5 * exp(-log(3) / 4 * y);
	// // vmax = 3.0 * exp(-log(3) / 6 * sqrt(x * x + y * y));

	// if (y >= 1.8 && y <= 2.0)
	// {
	// 	vmax *= 1.6;
	// }
	// else if (y > 2.0)
	// {
	// 	vmax *= 0.4;
	// }
	// else if (y > 1.4)
	// {
	// 	vmax *= 1.1;
	// }
	// // // vmax = 3 / 17 / 43 * (y * y - 60 * y) + 3;
	// // // vmax = 1.0 / 25.0 * (x - 7.5) * (x - 7.5) + 3.0 / 4.0;
	// // // vmax = var[3] * (x - 5.0) * (x - 5.0) / 25.0;
	// // // vmax = var[3];
	// // // vmax = 2.0 * (x - 5.0) * (x - 5.0) / 25.0 + 1.0;
	// // // vmax = 2.0 * abs(x-5.0)/5.0 + 1.0;
	// // // if (x >= 7.0 && x <= 8.0)
	// // // {
	// // // 	vx = vmax * 0.2;
	// // // }
	// // // else
	// // // vx = vmax * (1.0 - (y / 0.5) * (y / 0.5));
	// // // vx = vmax * (1.0 - (x / 0.5) * (x / 0.5));
	// vx = vmax;

	// // // if (x <= 4.7)
	// // // {
	// // // 	vx = var[3] * (-x * x / 4.7 / 4.7 * 0.2 + 0.8) * (1.0 - (y / 0.5) * (y / 0.5));
	// // // }
	// // // else if (x <= 5.3)
	// // // {
	// // // 	vx = var[3] * ((x - 4.7) * (x - 5.3) / 0.09 * 0.2 + 0.6) * (1.0 - (y / 0.5) * (y / 0.5));
	// // // }
	// // // else
	// // // {
	// // // 	vx = var[3] * (-(x - 10) * (x - 10) / 4.7 / 4.7 * 0.1 + 0.7) * (1.0 - (y / 0.5) * (y / 0.5));
	// // // }

	// result[0] = vx;

	// ! Diffusion desire function
	// result[0] = 64.0 * t * sin(PI * 2.0 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));

	// ! Convection desire function
	// * Function 1
	// if (x <= 0.5 && y <= 0.5)
	// {
	// 	// result[0] = pow(2 * x - 1.0, 2) * pow(2 * y - 1.0, 2);
	// 	result[0] = 1;
	// }
	// * Function 2

	// double theta = PI / 4 * 2.0 / 3.0;
	// double w_norm = 1.0 * 1.0;
	// vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};
	// double e = var[0];
	// result[0] = sin(t) * (x - (exp(w[0] * (x - 1) / e) - exp(-w[0] / e)) / (1 - exp(-w[0] / e))) * (y - (exp(w[1] * (y - 1) / e) - exp(-w[1] / e)) / (1 - exp(-w[1] / e)));
}

void UserSetting2D::SetVariables(string fn_par)
{
	var.resize(21);
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 21; i++)
		{
			fin >> stmp >> var[i];
		}
		fin.close();
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting2D::SetInitialCondition(string fn_in, string fn_out)
{
	int i, j;
	for (i = 0; i < state_num; i++)
		val_ini[i].resize(pts.size());

	if (ReadIC)
	{
		string fname(fn_in), stmp;
		stringstream ss;
		ifstream fin;
		fin.open(fname);
		if (fin.is_open())
		{
			for (j = 0; j < pts.size(); j++)
			{
				for (i = 0; i < state_num; i++)
				{
					fin >> stmp >> val_ini[i][j];
				}
			}
			fin.close();
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
		}
	}
	else
	{
		for (i = 0; i < pts.size(); i++)
		{
			val_ini[0][i] = 0.0;
		}
		TXTWriteIC(fn_in);
	}
	MPI_Barrier(comm);
	if (VisualizeIC && comRank == 0)
		VTKVisualizeIC(fn_out);
	PetscPrintf(PETSC_COMM_WORLD, "Set Initial Condition Done！\n");
}

void UserSetting2D::SetBoundaryCondition(string fn_in, string fn_out)
{
	int i, j;
	for (i = 0; i < 2; i++)
		val_bc[i].resize(pts.size());
	bc_flag.resize(pts.size(), 0);
	if (ReadBC)
	{
		string fname(fn_in), stmp;
		stringstream ss;
		ifstream fin;
		fin.open(fname);
		if (fin.is_open())
		{
			for (j = 0; j < pts.size(); j++)
			{
				for (i = 0; i < state_num; i++)
				{
					fin >> stmp >> val_bc[i][j];
				}
				fin >> stmp >> bc_flag[j];
			}
			fin.close();
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
		}
	}
	else
	{
		/// For one pipe
		double eps(1e-6);
		int count(0);
		// ! Diffusion test
		// for (i = 0; i < pts.size(); i++)
		// {

		// 	double x0 = pts[i].coor[0];
		// 	double x1 = pts[i].coor[1];
		// 	if (abs(x0 - 0.0) < eps || abs(x0 - 1.0) < eps || abs(x1 - 0.0) < eps || abs(x1 - 1.0) < eps)
		// 	{
		// 		if (x0 >= 0 & x0 <= 0.5 & x1 >= 0 & x1 <= 0.5)
		// 		{
		// 			// val_bc[0][i] = pow((2.0 * x0 - 1.0), 2.0) * pow((2.0 * x1 - 1.0), 2.0); // test1 and test2
		// 			val_bc[0][i] = 1.0; // test3 and test4
		// 		}
		// 		else
		// 		{
		// 			val_bc[0][i] = 0.0;
		// 		}
		// 		val_bc[0][i] = 0;
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}
		// 	// if (abs(x1 - 0.0) < eps || abs(x1 - 1.0) < eps)
		// 	// {
		// 	// 	bc_flag[i] = 0;
		// 	// }
		// 	// else if (abs(x0 - 0.0) < eps)
		// 	// {
		// 	// 	bc_flag[i] = 1;
		// 	// }
		// 	// else if (abs(x0 - 1.0) < eps)
		// 	// {
		// 	// 	bc_flag[i] = 2;
		// 	// }
		// 	// else
		// 	// {
		// 	// 	bc_flag[i] = -1;
		// 	// }

		// 	// if (abs(x0 - 0.0) < eps || abs(x0 - 1.0) < eps || abs(x1 - 0.0) < eps || abs(x1 - 1.0) < eps)
		// 	// {
		// 	// 	/// Concentration
		// 	// 	val_bc[0][i] = 0.0;
		// 	// 	bc_flag[i] = -1;
		// 	// 	continue;
		// 	// }
		// 	// if (abs(pts[i].coor[0] - 0.0) < eps)
		// 	// {
		// 	// 	/// Concentration
		// 	// 	val_bc[0][i] = 0.0;
		// 	// 	bc_flag[i] = -1;
		// 	// 	continue;
		// 	// }
		// 	// else if (abs(pts[i].coor[0] - 1.0) < eps)
		// 	// {
		// 	// 	/// Concentration
		// 	// 	val_bc[0][i] = 0.0;
		// 	// 	bc_flag[i] = -2;
		// 	// 	continue;
		// 	// }
		// 	// else if (abs(pts[i].coor[1] - 0.0) < eps)
		// 	// {
		// 	// 	val_bc[0][i] = 0.0;
		// 	// 	bc_flag[i] = -3;
		// 	// 	continue;
		// 	// }
		// 	// else if (abs(pts[i].coor[1] - 1.0) < eps)
		// 	// {
		// 	// 	val_bc[0][i] = 0.0;
		// 	// 	bc_flag[i] = -4;
		// 	// 	continue;
		// 	// }
		// 	// else if (abs(pts[i].coor[1] - 0.5) < eps || abs(pts[i].coor[1] + 0.5) < eps)
		// 	// {
		// 	// 	/// Concentration
		// 	// 	//	val_bc[0][i] = 1.0;			val_bc[1][i] = 0.0;			val_bc[2][i] = 2.0;
		// 	// 	/// v+
		// 	// 	val_bc[3][i] = 0.0;
		// 	// 	val_bc[4][i] = 0.0;
		// 	// 	/// v-
		// 	// 	val_bc[5][i] = 0.0;
		// 	// 	val_bc[6][i] = 0.0;
		// 	// 	bc_flag[i] = -3;
		// 	// 	continue;
		// 	// }
		// 	bc_flag[i] = count;

		// 	count++;
		// }
		// ! Convection test
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	if (abs(pts[i].coor[0] - 0.0) < eps || abs(pts[i].coor[0] - 1.0) < eps || abs(pts[i].coor[1] - 0.0) < eps || abs(pts[i].coor[1] - 1.0) < eps)
		// 	//if (abs(pts[i].coor[0] - 0.0) < eps || abs(pts[i].coor[0] - 1.0) < eps || abs(pts[i].coor[1] - 0.0) < eps )
		// 	{
		// 		bc_flag[i] = -1;
		// 		// if (abs(pts[i].coor[1] - 1.0) < eps || ((pts[i].coor[1] > 0.5) && (pts[i].coor[1] < 1.0) && abs(pts[i].coor[0] - 0.0) < eps))
		// 		// {
		// 		// 	/// Concentration
		// 		// 	val_bc[0][i] = 1.0;
		// 		// 	continue;
		// 		// }
		// 		// else
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	continue;
		// 		// }
		// 		// test 2
		// 		double x1 = pts[i].coor[0];
		// 		double x2 = pts[i].coor[1];
		// 		// if (x1 <= 0.5 && x2 <= 0.5)
		// 		// {
		// 		// 	val_bc[0][i] = pow(2 * x1 - 1.0, 2) * pow(2 * x2 - 1.0, 2);
		// 		// }
		// 		if (x1 > 0. && x1 <= 0.5 && abs(x2 - 0.0) < eps)
		// 		{
		// 			val_bc[0][i] = 1.0;
		// 		}
		// 		val_bc[0][i] = 0;
		// 		continue;
		// 	}
		// 	bc_flag[i] = count;
		// 	count++;
		// }
		// ! Single pipe model
		// eps = 1e-6;
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double vmax = 0.0, vx = 0.0, vy = 0.0;
		// 	vmax = var[3] * (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
		// 	vx = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));
		// 	// if (abs(pts[i].coor[0] - 10.0) < eps || abs(pts[i].coor[0] - 0.0) < eps )
		// 	// {
		// 	// 	val_bc[0][i] = 3.0;
		// 	// 	bc_flag[i] = -1;
		// 	// 	continue;
		// 	// }
		// 	if (abs(pts[i].coor[0] - 0.0) < eps)
		// 	{
		// 		// * simple test
		// 		val_bc[0][i] = 3.0;
		// 		// val_bc[0][i] = vx * 0.8;
		// 		// val_bc[0][i] = (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));

		// 		///// Concentration
		// 		//val_bc[0][i] = 1.0;
		// 		//val_bc[1][i] = 2.0;
		// 		//val_bc[2][i] = 0.0;
		// 		///// v+
		// 		//val_bc[3][i] = vx;
		// 		//val_bc[4][i] = vy;
		// 		//// val_bc[3][i] = 1.0;
		// 		//// val_bc[4][i] = 0.0;
		// 		///// v-
		// 		////val_bc[5][i] = -1.0;		val_bc[6][i] = 0.0;

		// 		///bc_flag
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}
		// 	else if (abs(pts[i].coor[0] - 10.0) < eps)
		// 	{
		// 		val_bc[0][i] = 1.0;
		// 		// val_bc[0][i] = vx * 0.7;
		// 		// val_bc[0][i] = (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));

		// 		///// Concentration
		// 		//val_bc[0][i] = 1.0;
		// 		//val_bc[1][i] = 2.0;
		// 		//val_bc[2][i] = 0.0;
		// 		///// v+
		// 		//val_bc[3][i] = vx;
		// 		//val_bc[4][i] = vy;
		// 		//// val_bc[3][i] = 1.0;
		// 		//// val_bc[4][i] = 0.0;
		// 		///// v-
		// 		////val_bc[5][i] = -1.0;		val_bc[6][i] = 0.0;

		// 		///bc_flag
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}
		// 	else if (abs(pts[i].coor[1] - 0.5) < eps || abs(pts[i].coor[1] + 0.5) < eps)
		// 	{
		// 		val_bc[0][i] = 0;
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}
		// 	// 	// else if (abs(pts[i].coor[0] - 10.0) < eps)
		// 	// 	// {
		// 	// 	// 	// * simple test
		// 	// 	// 	// val_bc[0][i] = 1.0;

		// 	// 	// 	val_bc[0][i] = vmax;
		// 	// 	// 	// /// Concentration
		// 	// 	// 	// val_bc[0][i] = 1.0;
		// 	// 	// 	// val_bc[1][i] = 0.0;
		// 	// 	// 	// val_bc[2][i] = 2.0;
		// 	// 	// 	// /// v+
		// 	// 	// 	// //val_bc[3][i] = 1.0;			val_bc[4][i] = 0.0;
		// 	// 	// 	// /// v-
		// 	// 	// 	// // val_bc[5][i] = -1.0;
		// 	// 	// 	// // val_bc[6][i] = 0.0;
		// 	// 	// 	// val_bc[5][i] = -vx;
		// 	// 	// 	// val_bc[6][i] = vy;
		// 	// 	// 	bc_flag[i] = -1; // ! Flag is changed
		// 	// 	// 	continue;
		// 	// 	// }
		// 	// 	// 	// if (abs(pts[i].coor[0] - 0.0) < eps || abs(pts[i].coor[0] - 10.0) < eps || abs(pts[i].coor[1] - 0.5) < eps || abs(pts[i].coor[1] + 0.5) < eps)
		// 	// 	// 	// {
		// 	// 	// 	// 	bc_flag[i] = -1;
		// 	// 	// 	// 	val_bc[0][i] = exp(-pts[i].coor[0] - 0.0 / var[0]) - 0.1 * pts[i].coor[0] - 0.0 + 2.0;
		// 	// 	// 	// 	continue;
		// 	// 	// 	// }
		// 	// 	// 	// else if (abs(pts[i].coor[1] - 0.5) < eps || abs(pts[i].coor[1] + 0.5) < eps)
		// 	// 	// 	// {
		// 	// 	// 	// 	/// Concentration
		// 	// 	// 	// 	//	val_bc[0][i] = 1.0;			val_bc[1][i] = 0.0;			val_bc[2][i] = 2.0;
		// 	// 	// 	// 	/// v+
		// 	// 	// 	// 	val_bc[3][i] = 0.0;
		// 	// 	// 	// 	val_bc[4][i] = 0.0;
		// 	// 	// 	// 	/// v-
		// 	// 	// 	// 	val_bc[5][i] = 0.0;
		// 	// 	// 	// 	val_bc[6][i] = 0.0;
		// 	// 	// 	// 	bc_flag[i] = -3;
		// 	// 	// 	// 	continue;
		// 	// 	// 	// }
		// 	bc_flag[i] = count;
		// 	count++;
		// }

		// ! Bifurcation model
		// eps = 0.048;
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double vmax = 0.0, vx = 0.0, vy = 0.0;
		// 	vmax = var[3];
		// 	vx = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));

		// 	double x = pts[i].coor[0], y = pts[i].coor[1];
		// 	double tt = sqrt(x * x + y * y);
		// 	double theta = atan2(y, x) * 180.0 / PI;
		// 	double tx, ty;

		// 	if (theta >= 120 || theta <= -120)
		// 	{
		// 		tx = -x;
		// 		ty = y;
		// 		// if(abs(abs(ty)-0.5)<eps)
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	bc_flag[i] = -1;
		// 		// 	continue;
		// 		// }
		// 		if (abs(tx - 3) < eps || abs(abs(ty) - 0.5) < eps)
		// 		{
		// 			// val_bc[0][i] = 3.0;
		// 			vmax = 3.0 * (tx + 3) / 6;
		// 			vmax = var[3];
		// 			val_bc[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 			bc_flag[i] = -1;
		// 			continue;
		// 		}
		// 	}
		// 	if (theta > 0 && theta < 120)
		// 	{
		// 		tx = x * 0.5 + y * 0.5 * sqrt(3);
		// 		ty = -x * 0.5 * sqrt(3) + y * 0.5;

		// 		if (abs(tx - 3) < eps || abs(abs(ty) - 0.5) < eps)
		// 		{
		// 			// val_bc[0][i] = 3.0;
		// 			vmax = var[3] * (-tx + 5) / 10;
		// 			// vmax = 3.0 * (-tx + 3)/6 + 1.0;
		// 			val_bc[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 			bc_flag[i] = -1;
		// 			continue;
		// 		}
		// 		// if(abs(abs(ty)-0.5)<eps)
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	bc_flag[i] = -1;
		// 		// 	continue;
		// 		// }
		// 		// if(abs(tx-3)<eps)
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	bc_flag[i] = -1;
		// 		// 	continue;
		// 		// }
		// 	}
		// 	if (theta <= 0 && theta > -120)
		// 	{
		// 		tx = x * 0.5 - y * 0.5 * sqrt(3);
		// 		ty = x * 0.5 * sqrt(3) + y * 0.5;

		// 		if (abs(tx - 3) < eps || abs(abs(ty) - 0.5) < eps)
		// 		{
		// 			// val_bc[0][i] = 3.0;
		// 			vmax = var[3] * (-tx + 5) / 10;
		// 			vmax = 3.0 * (-tx + 3) / 6;
		// 			vmax = 1.4 * var[3] * (-tx + 5) / 10;
		// 			// vmax = 1.0/6.0*(tx-1.0) *(tx-3.0) + 1.0;
		// 			val_bc[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 			bc_flag[i] = -1;
		// 			continue;
		// 		}
		// 		// if(abs(abs(ty)-0.5)<eps)
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	bc_flag[i] = -1;
		// 		// 	continue;
		// 		// }
		// 		// if(abs(tx-3)<eps)
		// 		// {
		// 		// 	val_bc[0][i] = 0.0;
		// 		// 	bc_flag[i] = -1;
		// 		// 	continue;
		// 		// }
		// 	}

		// 	if (abs(pts[i].coor[0] - (-3.0)) < eps)
		// 	{
		// 		// * simple test
		// 		// val_bc[0][i] = 3.0; // concentration
		// 		val_bc[0][i] = vmax; //velocity

		// 		///// Concentration
		// 		//val_bc[0][i] = 1.0;
		// 		//val_bc[1][i] = 2.0;
		// 		//val_bc[2][i] = 0.0;
		// 		///// v+
		// 		//val_bc[3][i] = vx;
		// 		//val_bc[4][i] = vy;
		// 		//// val_bc[3][i] = 1.0;
		// 		//// val_bc[4][i] = 0.0;
		// 		///// v-
		// 		////val_bc[5][i] = -1.0;		val_bc[6][i] = 0.0;

		// 		///bc_flag
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}

		// 	bc_flag[i] = count;
		// 	count++;
		// }

		//! Real Geometry single pipe
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	if (pts[i].label == 1)
		// 	{
		// 		// val_bc[0][i] = 3.0;
		// 		val_bc[0][i] = 3.0;
		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}
		// 	else if (pts[i].label == 0)
		// 	{
		// 		// val_bc[0][i] = 2.0;
		// 		val_bc[0][i] = 1.0;

		// 		// double y = pts[i].coor[1];
		// 		// val_bc[0][i] = 2.5 * exp(-log(3) / 4 * y);
		// 		// // vmax = 3.0 * exp(-log(3) / 6 * sqrt(x * x + y * y));

		// 		// if (y >= 1.8 && y <= 2.0)
		// 		// {
		// 		// 	val_bc[0][i] *= 1.6;
		// 		// }
		// 		// else if (y > 2.0)
		// 		// {
		// 		// 	val_bc[0][i] *= 0.4;
		// 		// }
		// 		// else if (y > 1.4)
		// 		// {
		// 		// 	val_bc[0][i] *= 1.1;
		// 		// }

		// 		bc_flag[i] = -1;
		// 		continue;
		// 	}

		// 	bc_flag[i] = count;
		// 	count++;
		// }

		//! Real Geometry
		for (i = 0; i < pts.size(); i++)
		{
			if (pts[i].label == 1)
			{
				val_bc[0][i] = 3.0;
				val_bc[0][i] = val_desire[0][i];
				bc_flag[i] = -1;
				continue;
			}
			// else if (pts[i].label == 0)
			// {
			// 	val_bc[0][i] = 0.0;
			// 	val_bc[0][i] = val_desire[0][i];
			// 	bc_flag[i] = -1;
			// 	continue;
			// }

			bc_flag[i] = count;
			count++;
		}

		TXTWriteBC(fn_in);
	}
	MPI_Barrier(comm);
	SetBoundaryMapping();

	if (VisualizeBC && comRank == 0)
		VTKVisualizeBC(fn_out);
	PetscPrintf(PETSC_COMM_WORLD, "Set Boundary Condition Done！\n");
}

void UserSetting2D::SetBoundaryMapping()
{
	int i, count(0);
	for (i = 0; i < bc_flag.size(); i++)
	{
		if (bc_flag[i] >= 0)
			nonbc_mapping.push_back(i);
		else
			bc_mapping.push_back(i);
	}
	n_bcpt = bc_mapping.size();
	PetscPrintf(PETSC_COMM_WORLD, "Number of non-BC pts: %d\n", nonbc_mapping.size());
	PetscPrintf(PETSC_COMM_WORLD, "Number of BC pts: %d\n", n_bcpt);

	// cout << "Number of non-BC pts: " << nonbc_mapping.size() << endl;
	// cout << "Number of BC pts: " << n_bcpt << endl;
}

void UserSetting2D::SetDesireState(string fn_in, string fn_out)
{
	int i, j;
	for (i = 0; i < 2; i++)
		val_desire[i].resize(pts.size());
	if (ReadDesire)
	{
		string fname(fn_in), stmp;
		stringstream ss;
		ifstream fin;
		fin.open(fname);
		if (fin.is_open())
		{
			for (j = 0; j < pts.size(); j++)
			{
				for (i = 0; i < state_num; i++)
				{
					fin >> val_desire[i][j];

					if (pts[j].label == -2)
						val_desire[i][j] = val_desire[i][j] * 5;
					if (pts[j].label == -3)
						val_desire[i][j] = val_desire[i][j] * 20;
					if (pts[j].label == -4)
						val_desire[i][j] = val_desire[i][j] * 4;
					else
						val_desire[i][j] = val_desire[i][j] * 1;
					// if (pts[j].coor[1] >= 1.8 && pts[j].coor[1] <= 2.0)
					// 	val_desire[i][j] = val_desire[i][j] * 0.3;
					// else if (pts[j].coor[1] > 2.0)
					// 	val_desire[i][j] = val_desire[i][j] * (0.3 + 0.4 / 1.7 * (pts[j].coor[1] - 2.0));
					// else if (pts[j].coor[1] > 1.2)
					// 	val_desire[i][j] = val_desire[i][j] * (0.4 + 0.6 / 0.8 * (2.0 - pts[j].coor[1]));

					// double x = pts[j].coor[0];
					// double y = pts[j].coor[1];
					// double z = pts[j].coor[2];
					// double d = -0.402 * x + 0.541 * y - 0.739 * z - 5.2573;
					// //velocity
					// // if (d < 0.54)
					// // 	val_desire[i][j] = val_desire[i][j] * 2.3;
					// // else if (d < 0.84 && j < 150)
					// // 	val_desire[i][j] = val_desire[i][j] * 2.5;
					// // if (j > 516 && j < 546)
					// // 	val_desire[i][j] = val_desire[i][j] * 0.3;
					// // else if ((j >= 547 && j < 1200) || j > 1232)
					// // 	val_desire[i][j] = val_desire[i][j] * 0.5;

					// // concentration
					// if (d < 0.54)
					// 	val_desire[i][j] = val_desire[i][j] * (2.6 - 1.5 / 12 * (sqrt(x * x + y * y + z * z) - 22)) / 1.2;
					// else if (d < 0.84 && j < 150)
					// 	val_desire[i][j] = val_desire[i][j] * 1.8;

					// if (j > 516 && j < 546)
					// 	val_desire[i][j] = val_desire[i][j] * 0.8;
					// else if ((j >= 547 && j < 1200) || j > 1232)
					// 	val_desire[i][j] = val_desire[i][j] * 0.8;
					// val_desire[i][j] *= 3.0;

					// if (j > 1343 && j < 1419)
					// 	val_desire[i][j] = val_desire[i][j] * 6.0;
				}
			}
			fin.close();
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
		}
	}
	else
	{
		/// For one pipe
		// double eps(1e-6);
		// double c1 = 0.40, c2 = 0.60;
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double x1 = pts[i].coor[0];
		// 	double x2 = pts[i].coor[1];

		// 	if ((pts[i].coor[0] >= c1 & pts[i].coor[0] <= c2) || (pts[i].coor[1] >= c1 & pts[i].coor[1] <= c2))
		// 	{
		// 		val_desire[0][i] = - pts[i].coor[0] * exp( - (pts[i].coor[0]-0.5)*(pts[i].coor[0]-0.5) - (pts[i].coor[1]-0.5)*(pts[i].coor[1]-0.5));
		// 	}
		// 	else
		// 	{
		// 		for(j = 0;j<state_num;j++)
		// 		{
		// 			val_desire[j][i] = 0.0;
		// 		}
		// 	}

		// 	if((pts[i].coor[0] >=  0 & pts[i].coor[0] <= 0.5) && (pts[i].coor[1] >=  0 & pts[i].coor[1] <= 0.5))
		// 	{
		// 		val_desire[0][i] = 1.0;
		// 	}
		// 	else
		// 	{
		// 		for(j = 0;j<state_num;j++)
		// 		{
		// 			val_desire[j][i] = 0.0;
		// 		}
		// 	}

		// 	//val_desire[0][i]= pts[i].coor[0]* (1.0 - pts[i].coor[0])* exp(-pts[i].coor[0]) * pts[i].coor[1]* (1.0 - pts[i].coor[1])* exp(-pts[i].coor[1]);
		// 	//  val_desire[0][i]=- pts[i].coor[0] *exp(-pow( pts[i].coor[0] - 0.5,2)-pow( pts[i].coor[1]-0.5,2));

		// 	for(j = 0;j<state_num;j++)
		// 	{
		// 		val_desire[j][i] = 0.0;
		// 	}
		// 	if(x1<=0.5 && x2<=0.5)
		// 	{
		// 		val_desire[0][i] = pow(2 * x1 - 1.0, 2) * pow(2 * x2 - 1.0, 2);
		// 	}
		// 	// if (x1>0 && x1 <= 0.5 )
		// 	// {
		// 	// 	val_desire[0][i] = 1.0;
		// 	// }
		// 	// for (j = 0; j < state_num; j++)
		// 	// {
		// 	// 	val_desire[j][i] = 0.0;
		// 	// }
		// }
		// //! For one pipe
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double vmax = 0.0;
		// 	double x0 = pts[i].coor[0];
		// 	// * simple test
		// 	// val_desire[0][i] = (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
		// 	val_desire[0][i] = 3.0 / 50.0 * x0 * x0 - 4.0 / 5.0 * x0 + 3.0;
		// 	//vmax = (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
		// 	//val_desire[3][i] = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));
		// 	//val_desire[4][i] = 0.0;
		// 	//val_desire[5][i] = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));
		// 	//val_desire[6][i] = 0.0;
		// }
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double x = pts[i].coor[0];
		// 	double y = pts[i].coor[1];
		// 	val_desire[0][i] = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0; // concentration
		// 																 // result = 3.0 / 50.0 * x0 * x0 - 4.0 / 5.0 * x0 + 1.0; // velocity
		// 																 // 	val_desire[0][i] = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0;
		// 																 // 	val_desire[0][i] = exp(-x / var[0]) - 0.1 * x + 2.0;
		// 	double vmax = 0.0, vx = 0.0, vy = 0.0;
		// 	// vmax = var[3] * (x - 5.0) * (x - 5.0) / 25.0;
		// 	vmax = var[3];
		// 	// vmax = var[3] * (10.2-x)/10.0;
		// 	if (x >= 7 && x <= 8)
		// 		vmax = 0.0;

		// 	vx = vmax * (1.0 - (y / 0.5) * (y / 0.5));

		// 	val_desire[0][i] = vx;
		// }
		//! For square domain
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double vmax = 0.0;
		// 	double x0 = pts[i].coor[0];
		// 	double x1 = pts[i].coor[1];
		// 	if (x0 >= 0 & x0 <= 0.5 & x1 >= 0 & x1 <= 0.5)
		// 	{
		// 		// val_desire[0][i] = pow((2.0 * x0 - 1.0), 2.0) * pow((2.0 * x1 - 1.0), 2.0); // test1 and test2
		// 		val_desire[0][i] = 1.0; // test3 and test4
		// 	}
		// 	else
		// 	{
		// 		val_desire[0][i] = 0.0;
		// 	}

		// 	// double theta = PI / 4 * 1.0 / 3.0;
		// 	// double w_norm = 1.0 * 2.0;
		// 	// vector<double> w = {w_norm * cos(theta), w_norm * sin(theta)};
		// 	// double e = var[0];
		// 	// val_desire[0][i] = sin(1) * (x0 - (exp(w[0] * (x0 - 1) / e) - exp(-w[0] / e)) / (1 - exp(-w[0] / e))) * (x1 - (exp(w[1] * (x1 - 1) / e) - exp(-w[1] / e)) / (1 - exp(-w[1] / e)));

		// 	// * simple test
		// 	// val_desire[0][i] = (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
		// 	// val_desire[0][i] = 64.0 * sin(PI * 2.0 * ((pts[i].coor[0] - 0.5) * (pts[i].coor[0] - 0.5) + (pts[i].coor[1] - 0.5) * (pts[i].coor[1] - 0.5)));
		// 	//vmax = (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
		// 	//val_desire[3][i] = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));
		// 	//val_desire[4][i] = 0.0;
		// 	//val_desire[5][i] = vmax * (1.0 - (pts[i].coor[1] / 0.5) * (pts[i].coor[1] / 0.5));
		// 	//val_desire[6][i] = 0.0;
		// }
		//! For bifurcation sym
		// for (i = 0; i < pts.size(); i++)
		// {
		// 	double x = pts[i].coor[0];
		// 	double y = pts[i].coor[1];
		// 	double tt = sqrt(x*x + y*y);
		// 	double tx, ty;
		// 	double vmax;
		// 	double theta = atan2(y,x) * 180.0 / PI;

		// 	val_desire[0][i] = 5.0;
		// 	if(theta >= 120 || theta <= -120)
		// 	{
		// 		tx = -x;
		// 		ty = y;

		// 		vmax = var[3];
		// 		val_desire[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 	}
		// 	if(theta > 0 && theta< 120)
		// 	{
		// 		tx = x * 0.5 + y * 0.5 * sqrt(3);
		// 		ty = -x * 0.5 * sqrt(3) + y * 0.5;

		// 		vmax = var[3];
		// 		val_desire[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 		// cout << "index: "<< i << " x: "<< x << " y: " << y << " desire: " << val_desire[0][i] << endl;
		// 	}
		// 	if(theta <= 0 && theta > -120)
		// 	{
		// 		tx = x * 0.5 - y * 0.5 * sqrt(3);
		// 		ty = x * 0.5 * sqrt(3) + y * 0.5;

		// 		vmax = var[3];
		// 		val_desire[0][i] = vmax * (1.0 - (ty / 0.5) * (ty / 0.5));
		// 	}
		// 	// if(theta >= 120 || theta <= -120)
		// 	// {
		// 	// 	val_desire[0][i] = var[3];
		// 	// }
		// 	// if(theta > 0 && theta< 120)
		// 	// {
		// 	// 	val_desire[0][i] =var[3] *( 1.0/6.0 * tt * tt -2.0/3.0 * tt + 1);
		// 	// }
		// 	// if(theta <= 0 && theta > -120)
		// 	// {
		// 	// 	val_desire[0][i] = var[3];
		// 	// }

		// 	// if(theta >= 120 || theta <= -120)
		// 	// {
		// 	// 	val_desire[0][i] = 2.0;
		// 	// }
		// 	// if(theta > 0 && theta< 120)
		// 	// {
		// 	// 	val_desire[0][i] =0.5;
		// 	// }
		// 	// if(theta <= 0 && theta > -120)
		// 	// {
		// 	// 	val_desire[0][i] = 1.0;
		// 	// }
		// }

		//! For real geometry pipe
		for (i = 0; i < pts.size(); i++)
		{
			double x = pts[i].coor[0];
			double y = pts[i].coor[1];
			double vmax = 0.0, vx = 0.0, vy = 0.0;
			// vmax = 3.0 / 50.0 * x * x - 4.0 / 5.0 * x + 3.0;
			// vmax = 1.0 / 70.0 * (x - 12) * (x - 12) + 33 / 35 + 0.2;
			// // vmax = 1.0 + (10 - x) / 5;
			vmax = 2.5 * exp(-log(3) / 4 * y);
			// vmax = 3.0 * exp(-log(3) / 6 * sqrt(x * x + y * y));

			if (y >= 1.8 && y <= 2.0)
			{
				vmax *= 1.6;
			}
			else if (y > 2.0)
			{
				vmax *= 0.4;
			}
			else if (y > 1.4)
			{
				vmax *= 1.1;
			}
			vx = vmax;

			val_desire[0][i] = vx;
		}

		TXTWriteDesire(fn_in);
	}
	MPI_Barrier(comm);
	if (VisualizeDesire && comRank == 0)
		VTKVisualizeDesire(fn_out);
	PetscPrintf(PETSC_COMM_WORLD, "Set Boundary Condition Done！\n");
}

void UserSetting2D::TXTWriteIC(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		for (j = 0; j < pts.size(); j++)
		{
			for (i = 0; i < state_num; i++)
			{

				if (i == 1)
				{
					fout << val_ini[i][j] << "\n";
				}
				else
				{
					fout << val_ini[i][j] << " ";
				}
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::TXTWriteBC(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		for (j = 0; j < pts.size(); j++)
		{
			for (i = 0; i < state_num; i++)
			{
				fout << val_bc[i][j] << " ";
				// if (i == 6)
				// {
				// 	fout << val_bc[i][j] << "\n";
				// }
				// else
				// {
				// 	fout << val_bc[i][j] << " ";
				// }
			}
			fout << bc_flag[j] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::TXTWriteDesire(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		for (j = 0; j < pts.size(); j++)
		{
			for (i = 0; i < state_num; i++)
			{
				if (i != 1)
					fout << val_desire[i][j] << " ";
				else
					fout << val_desire[i][j] << "\n";
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::VTKVisualizeIC(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (i = 0; i < pts.size(); i++)
		{
			fout << std::setprecision(9) << pts[i].coor[0] << std::fixed << " " << std::setprecision(9) << pts[i].coor[1] << std::fixed << " " << std::setprecision(9) << pts[i].coor[2] << std::fixed << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << "\n";
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\n";
		fout << "VECTORS y float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[0][i] << " " << 0 << " " << 0 << "\n";
		// fout << "VECTORS u float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[2][i] << " " << val_ini[3][i] << " " << 0 << "\n";
		// fout << "VECTORS lambda float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[4][i] << " " << val_ini[5][i] << " " << 0 << "\n";
		// fout << "SCALARS N0 float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[0][i] << "\n";
		// fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[1][i] << "\n";
		// fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[2][i] << "\n";
		// fout << "VECTORS Vplus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[3][i] << " " << val_ini[4][i] << " " << 0 << "\n";
		// fout << "VECTORS Vminus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[5][i] << " " << val_ini[6][i] << " " << 0 << "\n";
		// fout << "VECTORS Fplus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[7][i] << " " << val_ini[8][i] << " " << 0 << "\n";
		// fout << "VECTORS Fminus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_ini[9][i] << " " << val_ini[10][i] << " " << 0 << "\n";
		// for (j = 0; j < 7; j++)
		// {
		// 	fout << "SCALARS lambda" << j+1 << " float 1\nLOOKUP_TABLE default\n";
		// 	for (i = 0; i < pts.size(); i++)
		// 		fout << val_ini[j + 11][i] << "\n";
		// }

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::VTKVisualizeBC(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (i = 0; i < pts.size(); i++)
		{
			fout << std::setprecision(9) << pts[i].coor[0] << std::fixed << " " << std::setprecision(9) << pts[i].coor[1] << std::fixed << " " << std::setprecision(9) << pts[i].coor[2] << std::fixed << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\n";
		fout << "SCALARS bc_flag float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << bc_flag[i] << "\n";
		// fout << "SCALARS N0 float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_bc[0][i] << "\n";
		// fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_bc[1][i] << "\n";
		// fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_bc[2][i] << "\n";
		fout << "VECTORS y float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[0][i] << " " << 0 << " " << 0 << "\n";
		// fout << "VECTORS Vminus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_bc[5][i] << " " << val_bc[6][i] << " " << 0 << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::VTKVisualizeDesire(string fn_out)
{
	string fname = fn_out;
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (i = 0; i < pts.size(); i++)
		{
			fout << std::setprecision(9) << pts[i].coor[0] << std::fixed << " " << std::setprecision(9) << pts[i].coor[1] << std::fixed << " " << std::setprecision(9) << pts[i].coor[2] << std::fixed << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i < mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\n";
		// fout << "SCALARS N0 float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_desire[0][i] << "\n";
		// fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_desire[1][i] << "\n";
		// fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_desire[2][i] << "\n";
		fout << "VECTORS yd float\n";
		for (i = 0; i < pts.size(); i++)
			fout << 0 << " " << 0 << " " << val_desire[0][i] << "\n";
		// fout << "VECTORS Vminus float\n";
		// for (i = 0; i < pts.size(); i++)
		// 	fout << val_desire[5][i] << " " << val_desire[6][i] << " " << 0 << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::ReadMesh(string fn)
{
	string fname(fn), stmp, stmp_line;
	int npts, neles, itmp;
	double dtmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		while (getline(fin, stmp))
		{
			if (fin)
			{
				string upper(stmp);

				std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
				if (upper.find("POINTS") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> npts >> stmp_line;
					nPoint = npts;
					PetscPrintf(PETSC_COMM_WORLD, "Start reading node # of points: %d\n", npts);
					pts.resize(npts);
					for (int i = 0; i < npts; i++)
					{
						fin >> pts[i].coor[0] >> pts[i].coor[1] >> pts[i].coor[2];
					}
					//cout << "Node number: " << pts.size() << endl;
				}
				else if (upper.find("CELLS") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> neles >> stmp_line;
					PetscPrintf(PETSC_COMM_WORLD, "Start reading cells # of elements: %d\n", neles);
					mesh.resize(neles);
					for (int i = 0; i < neles; i++)
					{
						fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
						for (int j = 0; j < 4; j++)
						{
							mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
							mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
							mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
						}
					}
					//cout << "Node number: " << pts.size() << endl;
				}
				else if (upper.find("POINT_DATA") != string::npos)
				{
					stringstream sstmp(upper);
					sstmp >> stmp_line >> npts;
					// cout << "Start reading point labels # of labels: " << npts << endl;
					PetscPrintf(PETSC_COMM_WORLD, "Start reading point labels # of labels: %d\n", npts);
					for (int i = 0; i < 2; i++)
					{
						getline(fin, stmp); //skip lines
					}
					for (int i = 0; i < npts; i++)
					{
						fin >> pts[i].label;
					}
					//cout << "Node number: " << cp.size() << endl;
				}
				// 2.) Look for the "*Element" section
				continue;
			}
		}
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}

	// if (fin.is_open())
	// {
	// 	for (int i = 0; i < 4; i++)
	// 		getline(fin, stmp); //skip lines
	// 	fin >> stmp >> npts >> stmp;
	// 	pts.resize(npts);
	// 	for (int i = 0; i < npts; i++)
	// 	{
	// 		fin >> pts[i].coor[0] >> pts[i].coor[1] >> pts[i].coor[2];
	// 	}
	// 	getline(fin, stmp);
	// 	fin >> stmp >> neles >> itmp;
	// 	mesh.resize(neles);
	// 	for (int i = 0; i < neles; i++)
	// 	{
	// 		fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
	// 		for (int j = 0; j < 4; j++)
	// 		{
	// 			mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
	// 			mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
	// 			mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
	// 		}
	// 	}
	// 	// 		ofstream fout;
	// 	// fout.open(work_dir + "debug/readmesh.txt");
	// 	for (int i = 0; i < neles + 6; i++)
	// 	{
	// 		getline(fin, stmp); //skip lines

	// 		// fout << stmp <<endl;
	// 	}

	// 	for (int i = 0; i < npts; i++)
	// 	{
	// 		fin >> pts[i].label;
	// 		// fout << pts[i].label <<endl;
	// 	}

	// 	fin.close();
	// 	PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	// 	// PetscPrintf(PETSC_COMM_WORLD, "# of Nodes: %d \n", pts.size());
	// 	// PetscPrintf(PETSC_COMM_WORLD, "# of elements: %d \n", mesh.size());
	// }
	// else
	// {
	// 	PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	// }
}

// void UserSetting2D::CreateDMPlexMesh(string fn)
// {
// 	vector<vector<int>> bzmesh_dmplex;
// 	string fname(fn), stmp;

// 	int npts, neles, itmp;
// 	double dtmp;
// 	ifstream fin;
// 	fin.open(fname);
// 	if (fin.is_open())
// 	{
// 		fin >> neles;
// 		bzmesh_dmplex.resize(neles);
// 		getline(fin, stmp);
// 		for (int i = 0; i < neles; i++)
// 		{
// 			getline(fin, stmp);
// 			stringstream sstream(stmp);
// 			while (sstream >> itmp)
// 			{
// 				//cout << itmp << endl;
// 				bzmesh_dmplex[i].push_back(itmp - 1 + neles);
// 			}
// 		}
// 		fin.close();
// 		PetscPrintf(PETSC_COMM_WORLD, "Bezier Mesh For DMPlex Loaded!\n");
// 		// PetscPrintf(PETSC_COMM_WORLD, "# of Nodes: %d \n", pts.size());
// 		// PetscPrintf(PETSC_COMM_WORLD, "# of elements: %d \n", mesh.size());

// 		// * Create DMPlex togology
// 		DMPlexCreate(PETSC_COMM_WORLD, &dm_bzmesh);
// 		DMPlexSetChart(dm_bzmesh, 0, neles + pts.size());
// 		for (int i = 0; i < neles; i++)
// 		{
// 			DMPlexSetConeSize(dm_bzmesh, i, bzmesh_dmplex[i].size());
// 		}

// 		DMSetUp(dm_bzmesh);
// 		for (int i = 0; i < neles; i++)
// 		{
// 			const PetscInt *dm_IEN = &bzmesh_dmplex[i][0];
// 			DMPlexSetCone(dm_bzmesh, i, dm_IEN);
// 		}
// 		DMSetDimension(dm_bzmesh, 1);
// 		DMPlexSymmetrize(dm_bzmesh);
// 		DMPlexStratify(dm_bzmesh);

// 		for (int i = 0; i < neles; i++)
// 		{
// 			DMPlexSetCellType(dm_bzmesh, i, DM_POLYTOPE_SEGMENT);
// 		}

// 		DMViewFromOptions(dm_bzmesh, NULL, "-dm_view");
// 		PetscPrintf(PETSC_COMM_WORLD, "DMPlex mesh created!\n");

// 		// * Set up DMPlex PetscSection
// 		PetscInt numFields = 1;
// 		PetscInt numDOF[numFields * (1 + 1)], numComp[numFields];
// 		PetscInt numBC = 1;
// 		PetscInt bcField[numBC];
// 		PetscSection section;

// 		/*   Number of Components per Field  */
// 		nTstep = 1;
// 		for (int i = 0; i < numFields; i++)
// 		{
// 			numComp[i] = 1 * nTstep;
// 		}

// 		/*   Spatial DOF of field component (0 to not use)   */
// 		for (int f = 0; f < numFields; f++)
// 		{
// 			numDOF[f * (1 + 1) + 0] = 1;
// 			numDOF[f * (1 + 1) + 1] = 0;
// 		}
// 		// for (int d = 0; d < 1 + 1; d++)
// 		// 	numDOF[f * (1 + 1) + d] = 0;

// 		/*   Specificy DOF of field component that we intend to use  */
// 		// numDOF[0] = 1;
// 		// numDOF[2] = 1;
// 		// numDOF[6] = 1;

// 		/*   Boundary condition of field component set to Dirichlet  */
// 		// bcField[0] = 0;

// 		ierr = DMSetNumFields(dm_bzmesh, numFields);
// 		ierr = DMPlexCreateSection(dm_bzmesh, NULL, numComp, numDOF, NULL, NULL, NULL, NULL, NULL, &section);
// 		ierr = PetscSectionSetFieldName(section, 0, "State");
// 		// ierr = PetscSectionSetFieldName(section, 1, "Ctrl");
// 		// ierr = PetscSectionSetFieldName(section, 2, "Lambda");
// 		ierr = DMSetSection(dm_bzmesh, section);

// 		/*   Dont forget to clean up */
// 		ierr = PetscSectionDestroy(&section);
// 		PetscPrintf(PETSC_COMM_WORLD, "DMPlex section created!\n");

// 		Mat mat_test;
// 		DMSetMatType(dm_bzmesh, MATAIJ);
// 		DMCreateMatrix(dm_bzmesh, &mat_test);
// 		DebugVisualizeMat(mat_test, work_dir + "debug/dm_mat.m");

// 		PetscErrorCode ierr;
// 		PetscViewer viewer;
// 		PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

// 		string fname_dm = work_dir + "debug/dm_bzmesh.m";

// 		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname_dm.c_str(), &viewer);
// 		ierr = PetscViewerPushFormat(viewer, format);
// 		ierr = DMView(dm_bzmesh, viewer);
// 		ierr = PetscViewerPopFormat(viewer);
// 		ierr = PetscViewerDestroy(&viewer);
// 	}
// 	else
// 	{
// 		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
// 	}
// }

void UserSetting2D::AssignProcessor(string fn)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while (!fin.eof() && fin.peek() != EOF)
		{
			fin >> tmp;
			ele_process[tmp].push_back(i);
			i++;
			fin.get();
		}
		n_bzmesh = i;
		PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished!\n");
		fin.close();
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting2D::Readbzmeshinfo(string fn)
{
	string fname(fn), stmp;
	int neles, itmp;
	vector<int> IEN_tmp;

	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		fin >> neles;
		bzmeshinfo.resize(nPoint);
		getline(fin, stmp);
		for (int i = 0; i < neles; i++)
		{
			IEN_tmp.clear();
			getline(fin, stmp);
			stringstream sstream(stmp);
			while (sstream >> itmp)
			{
				IEN_tmp.push_back(itmp - 1);
			}
			for (int j = 0; j < IEN_tmp.size(); j++)
			{
				bzmeshinfo[IEN_tmp[j]].insert(bzmeshinfo[IEN_tmp[j]].end(), IEN_tmp.begin(), IEN_tmp.end());
			}
		}
		fin.close();
		for (int i = 0; i < bzmeshinfo.size(); i++)
		{
			unordered_set<int> s;
			for (int j : bzmeshinfo[i])
				s.insert(j);
			bzmeshinfo[i].assign(s.begin(), s.end());
			sort(bzmeshinfo[i].begin(), bzmeshinfo[i].end());
		}

		ofstream fout;
		fout.open(work_dir + "debug/bzmesh_cnct.txt");
		for (auto vec : bzmeshinfo)
		{
			for (auto i : vec)
			{
				fout << i << " ";
			}
			fout << "\n";
		}

		PetscPrintf(PETSC_COMM_WORLD, "Bezier Mesh information Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}
