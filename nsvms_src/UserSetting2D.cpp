#include "UserSetting2D.h"

UserSetting2D::UserSetting2D()
{
}

void UserSetting2D::SetVariables(string fn_par, vector<double> &var)
{
	var.resize(12);
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 12; i++)
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

void UserSetting2D::SetInitialCondition(int ndof, vector<double> &Vel0, vector<double> &Pre0, vector<Vertex2D> &pts, const vector<array<double, dim>> velocity_node, const int dim)
{
	double val_ini[2] = {0.0, 0.0}; //Vel0, Pre0
	Vel0.clear();
	Pre0.clear();
	Vel0.resize(ndof * dim, val_ini[0]);
	Pre0.resize(ndof, val_ini[1]);
}

void UserSetting2D::ReadMesh(string fn, vector<Vertex2D> &pts, vector<Element2D> &mesh)
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
					// cout << "Start reading node # of points: " << npts << endl;
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
					// cout << "Start reading cells # of elements: " << neles << endl;
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
					for (int i = 0; i < 2; i++)
					{
						getline(fin, stmp); //skip lines
					}
					for (int i = 0; i < npts; i++)
					{
						fin >> pts[i].label;
					}
					//cout << "Node number: " << pts.size() << endl;
				}
				// 2.) Look for the "*Element" section
				continue;
			}
			//if (fin.eof())
			//{
			//	cout << "Read to the end!" << endl;
			//	break;
			//}
		}
		fin.close();
		PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
		// for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		// fin >> stmp >> npts >> stmp;
		// pts.resize(npts);
		// for (int i = 0; i < npts; i++)
		// {
		// 	fin >> pts[i].coor[0] >> pts[i].coor[1] >> pts[i].coor[2];
		// }
		// getline(fin, stmp);
		// fin >> stmp >> neles >> itmp;
		// mesh.resize(neles);
		// for (int i = 0; i < neles; i++)
		// {
		// 	fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
		// 	for (int j = 0; j < 4; j++)
		// 	{
		// 		mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
		// 		mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
		// 	}

		// }
		// for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		// for (int i = 0; i < npts; i++)	fin >> pts[i].label;
		// fin.close();
		// PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting2D::ReadVelocityField(string fn, int npts, vector<array<double, dim>> &velocity)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	velocity.resize(npts);
	double dtmp;
	if (fin.is_open())
	{
		for (int i = 0; i < npts; i++)
		{
			fin >> velocity[i][0] >> velocity[i][1] >> dtmp;
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