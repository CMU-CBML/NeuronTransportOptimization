#include "UserSetting2D.h"

UserSetting2D::UserSetting2D()
{
}

void UserSetting2D::InitializeUserSetting(string path)
{
	int rank, nProcs;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

	ele_process.resize(nProcs);
	
	work_dir = path;
	string fn_mesh(path + "controlmesh.vtk");
	string fn_bz(path + "bzmeshinfo.txt.epart." + to_string(nProcs));
	string fn_IC(path + "InitialCondition.txt");
	string fn_IC_visual(path + "InitialCondition.vtk");
	string fn_BC(path + "BoundaryCondition.txt");
	string fn_BC_visual(path + "BoundaryCondition.vtk");
	string fn_Desire(path + "DesireState.txt");
	string fn_Desire_visual(path + "DesireState.vtk");

	string fn_velocity(path + "initial_velocityfield.txt");		
	string fn_parameter(path + "simulation_parameter.txt");
	
	SetVariables(fn_parameter);

	// SetInitialCondition(int ndof, vector<double> &Vel0, vector<double> &Pre0, vector<Vertex2D> &pts, const vector<array<double, 2>> velocity_node);
	ReadMesh(fn_mesh);
	SetInitialCondition(fn_IC, fn_IC_visual);
	SetBoundaryCondition(fn_BC, fn_BC_visual);
	SetDesireState(fn_Desire,fn_Desire_visual);
	// ReadVelocityField(fn_velocity);
	AssignProcessor(fn_bz);

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
	for (i = 0; i < 6; i++)
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
				for (i = 0; i < 6; i++)
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
		for(i = 0; i < pts.size();i++)
		{
			/// Concentration
			val_ini[0][i] = 1.0;			val_ini[1][i] = 2.0;			val_ini[2][i] = 2.0;
			/// v+
			val_ini[3][i] = 1.0;			val_ini[4][i] = 0.0;
			/// v-
			val_ini[5][i] = -1.0;			val_ini[6][i] = 0.0;
			/// f+, f-
			val_ini[7][i] = 1.0;			val_ini[8][i] = 0.0;
			val_ini[9][i] = 1.0;			val_ini[10][i] = 0.0;
			/// lambda
			for(j=0;j<7;j++)
				val_ini[j+11][i] = 1.0;
		}
		TXTWriteIC(fn_in);
	}
	if(VisualizeIC)
		VTKVisualizeIC(fn_out);
	PetscPrintf(PETSC_COMM_WORLD, "Set Initial Condition Done！\n");
}

void UserSetting2D::SetBoundaryCondition(string fn_in,  string fn_out)
{
	int i, j;
	for (i = 0; i < 7; i++)
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
				for (i = 0; i < 7; i++)
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
		for (i = 0; i < pts.size(); i++)
		{
			double vmax = 0.0, vx = 0.0, vy =0.0;
			vmax = var[1] * (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
			vx = vmax * (1.0 - (pts[i].coor[1]/0.5) * (pts[i].coor[1]/0.5));
			if (abs(pts[i].coor[0] - 0.0) < eps)
			{
				/// Concentration
				val_bc[0][i] = 1.0;
				val_bc[1][i] = 2.0;
				val_bc[2][i] = 0.0;
				/// v+
				val_bc[3][i] = vx;
				val_bc[4][i] = vy;
				// val_bc[3][i] = 1.0;
				// val_bc[4][i] = 0.0;
				/// v-
				//val_bc[5][i] = -1.0;		val_bc[6][i] = 0.0;
				///bc_flag
				bc_flag[i] = -1;
				continue;
			}
			else if (abs(pts[i].coor[0] - 10.0) < eps)
			{
				/// Concentration
				val_bc[0][i] = 1.0;
				val_bc[1][i] = 0.0;
				val_bc[2][i] = 2.0;
				/// v+
				//val_bc[3][i] = 1.0;			val_bc[4][i] = 0.0;
				/// v-
				// val_bc[5][i] = -1.0;
				// val_bc[6][i] = 0.0;
				val_bc[5][i] = -vx;
				val_bc[6][i] = vy;
				bc_flag[i] = -2;
				continue;
			}
			else if (abs(pts[i].coor[1] - 0.5) < eps || abs(pts[i].coor[1] + 0.5) < eps)
			{
				/// Concentration
				//	val_bc[0][i] = 1.0;			val_bc[1][i] = 0.0;			val_bc[2][i] = 2.0;
				/// v+
				val_bc[3][i] = 0.0;
				val_bc[4][i] = 0.0;
				/// v-
				val_bc[5][i] = 0.0;
				val_bc[6][i] = 0.0;
				bc_flag[i] = -3;
				continue;
			}
			bc_flag[i] = count;
			count++;
		}
		TXTWriteBC(fn_in);
	}

	SetBoundaryMapping();

	if(VisualizeBC)
		VTKVisualizeBC(fn_out);
	PetscPrintf(PETSC_COMM_WORLD, "Set Boundary Condition Done！\n");
}

void UserSetting2D::SetBoundaryMapping()
{
	int i, count(0);
	for(i=0;i<bc_flag.size();i++)
	{
		if(bc_flag[i] >= 0)
			nonbc_mapping.push_back(i);
		else
			bc_mapping.push_back(i);
	}
	n_bcpt = bc_mapping.size();
	cout << "Number of non-BC pts: " << nonbc_mapping.size()<< endl;
	cout << "Number of BC pts: " << n_bcpt<< endl;
}

void UserSetting2D::SetDesireState(string fn_in,  string fn_out)
{
	int i, j;
	for (i = 0; i < 7; i++)
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
				for (i = 0; i < 7; i++)
				{
					fin >> stmp >> val_desire[i][j];
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
		double eps(1e-6);
		for (i = 0; i < pts.size(); i++)
		{
			double vmax = 0.0;
			vmax = (pts[i].coor[0] - 5.0) * (pts[i].coor[0] - 5.0) / 25.0;
			val_desire[3][i] = vmax * (1.0 - (pts[i].coor[1]/0.5) * (pts[i].coor[1]/0.5));
			val_desire[4][i] = 0.0;
			val_desire[5][i] = vmax * (1.0 - (pts[i].coor[1]/0.5) * (pts[i].coor[1]/0.5));
			val_desire[6][i] = 0.0;			
		}
		TXTWriteDesire(fn_in);
	}
	if(VisualizeDesire)
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
			for (i = 0; i < 18; i++)
			{

				if (i == 17)
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
			for (i = 0; i < 7; i++)
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
			for (i = 0; i < 7; i++)
			{
				if(i!=6)
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
		for (i = 0; i<pts.size(); i++)
		{
			fout << pts[i].coor[0] << " " << pts[i].coor[1] << " " << pts[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\nSCALARS N0 float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[0][i] << "\n";
		fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[1][i] << "\n";
		fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[2][i] << "\n";
		fout << "VECTORS Vplus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[3][i] << " " << val_ini[4][i] << " " << 0 << "\n";
		fout << "VECTORS Vminus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[5][i] << " " << val_ini[6][i] << " " << 0 << "\n";
		fout << "VECTORS Fplus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[7][i] << " " << val_ini[8][i] << " " << 0 << "\n";
		fout << "VECTORS Fminus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_ini[9][i] << " " << val_ini[10][i] << " " << 0 << "\n";
		for (j = 0; j < 7; j++)
		{
			fout << "SCALARS lambda" << j+1 << " float 1\nLOOKUP_TABLE default\n";
			for (i = 0; i < pts.size(); i++)
				fout << val_ini[j + 11][i] << "\n";
		}
		
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
		for (i = 0; i<pts.size(); i++)
		{
			fout << pts[i].coor[0] << " " << pts[i].coor[1] << " " << pts[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\n";
		fout << "SCALARS bc_flag float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << bc_flag[i] << "\n";
		fout << "SCALARS N0 float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[0][i] << "\n";
		fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[1][i] << "\n";
		fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[2][i] << "\n";
		fout << "VECTORS Vplus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[3][i] << " " << val_bc[4][i] << " " << 0 << "\n";
		fout << "VECTORS Vminus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_bc[5][i] << " " << val_bc[6][i] << " " << 0 << "\n";
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
		for (i = 0; i<pts.size(); i++)
		{
			fout << pts[i].coor[0] << " " << pts[i].coor[1] << " " << pts[i].coor[2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << pts.size() << "\n";
		fout << "SCALARS N0 float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_desire[0][i] << "\n";
		fout << "SCALARS Nplus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_desire[1][i] << "\n";
		fout << "SCALARS Nminus float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_desire[2][i] << "\n";
		fout << "VECTORS Vplus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_desire[3][i] << " " << val_desire[4][i] << " " << 0 << "\n";
		fout << "VECTORS Vminus float\n";
		for (i = 0; i < pts.size(); i++)
			fout << val_desire[5][i] << " " << val_desire[6][i] << " " << 0 << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void UserSetting2D::ReadMesh(string fn)
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	double dtmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> dtmp;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
			for (int j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
			}

		}
		for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		for (int i = 0; i < npts; i++)	fin >> pts[i].label;
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting2D::ReadVelocityField(string fn)
{
	
}

void UserSetting2D::AssignProcessor(string fn)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while(!fin.eof() && fin.peek() != EOF)
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


void UserSetting2D::SetVariables(string fn_par, vector<double>& var)
{
	var.resize(19);
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 19; i++)
		{
			fin >> stmp >> var[i];
		}	
		fin.close();
	}	
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
	//double var0[2] = { 1.0,0.5 };	//nu, rou
	//var.clear();
	//var.resize(2);
	//var.assign(var0, var0 + 2);
}

void UserSetting2D::SetInitialCondition(int ndof, vector<double>& Vel0, vector<double>& Pre0, vector<Vertex2D>& pts, const vector<array<double,2>> velocity_node)
{
	double val_ini[2] = { 0.0,0.0 };//Vel0, Pre0
	Vel0.clear();
	Pre0.clear();
	Vel0.resize(ndof * 2, val_ini[0]);
	Pre0.resize(ndof, val_ini[1]);
}

void UserSetting2D::ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh)
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	double dtmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> dtmp;
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
			for (int j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
			}

		}
		for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		for (int i = 0; i < npts; i++)	fin >> pts[i].label;
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}

}

void UserSetting2D::ReadVelocityField(string fn, int npts, vector<array<double,2>>& velocity)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	velocity.resize(npts);
	if (fin.is_open())
	{
		for (int i = 0; i < npts; i++)
		{
			for(int j = 0; j < 2; j++)
			{
				fin >> velocity[i][0] >> velocity[i][1] >> velocity[i][2];
			}
			
		}
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Velocity Field Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting2D::AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while(!fin.eof() && fin.peek() != EOF)
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