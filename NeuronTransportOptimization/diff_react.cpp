#include "diff_react.h"

void SetVariables(vector<double>& var)
{
	double var0[7] = { 0.1,5.0,-0., //Dn0, v_plus, v_minus
						1.0,0.0,	// k+, k-
						0.5,0.0	};	//k'+,k'-
	var.clear();
	var.resize(7);
	var.assign(var0, var0 + 7);
}

void SetInitialCondition(int ndof, int ndof_b, vector<double>& CA0, vector<double>& NX0, vector<double>& NB0, vector<array<double, 3>>& pts, const vector<int>& label, const vector<int>& pid_loc)
{
	double val_ini[3] = { 0.0,0.0,0.0 };//CA, N_plus, N_minus
	CA0.clear();
	NX0.clear();
	NB0.clear();
	CA0.resize(ndof, val_ini[0]);
	NX0.resize(ndof, val_ini[1]);
	NB0.resize(ndof, val_ini[2]);

	for (int i = 0; i < ndof; i++)
	{
		
		////cylinder ICs
		//if (label[i] == 1 && (sqrt(pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]) < 0.49))//income
		//{
		//	CA0[i] = 1.0;
		//	NX0[i] = 2.0;
		//}
		
		//cylinder photoactivation ICs
		if (pts[i][0]>4 && pts[i][0]<6)// && (sqrt(pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]) < 0.49))//income
		{
			CA0[i] = 1.0;
			NX0[i] = 2.0;
		}

		////bifurcation photoactivation ICs
		//if (i>602 && i<1005)// && (sqrt(pts[i][1] * pts[i][1] + pts[i][2] * pts[i][2]) < 0.49))//income
		//{
		//	CA0[i] = 1.0;
		//	NX0[i] = 2.0;
		//}

		//bifurcation ICs
		//double R0 = sqrt(pow(pts[i][0] - 0.0, 2) + pow(pts[i][1] - 0.0, 2) + pow(pts[i][2] - 0.0, 2));
		//if (label[i] == 1)
		//{
		//	CA0[i] = 1.0;
		//	NX0[i] = 2.0;
		//
		//}
		//if (label[i] > 1)
		//{
		//	NB0[i] = 0.0;
		//}

		//3 bifurcation ICs
		//double R0 = sqrt(pow(pts[i][0] - 0.0, 2) + pow(pts[i][1] - 0.0, 2) + pow(pts[i][2] - 0.0, 2));
		//if (label[i] == 1)
		//{
		//	if (R0 > 2.0)//inlet wall
		//	{
		//		CA0[i] = 0.0;
		//		NX0[i] = 0.0;
		//	}
		//	else//inlet
		//	{
		//		CA0[i] = 1.0;
		//		NX0[i] = 2.0;
		//	}
		//
		//}
		//if (label[i] > 1)
		//{
		//	NB0[i] = 0.0;
		//}
	}
}

void SetTempData(int ndof, int ndof_b, vector<double>& CA, vector<double>& NX, vector<double>& NB)
{
	CA.clear();
	NX.clear();
	NB.clear();
	CA.resize(ndof);
	NX.resize(ndof);
	NB.resize(ndof);
}

void SetTestBC(double& CAi, double& CAs, double& length)
{
	double val_bc[3] = { 1.0,0.3,1.};
	CAi = val_bc[0];
	CAs = val_bc[1];
	length = val_bc[2];
}

double MaxDifference(vector<double> a, vector<double> b)
{
	double c(0.0);
	vector<double> diff(a.size());
	for (int i = 0; i < a.size(); i++)
	{
		diff[i] = abs(a[i] - b[i])/abs(a[i]);
		if (diff[i] >= c)
		{
			c = diff[i];
		}
	}
	return c;
}

void ReadMesh(string fn, vector<array<double, 3>>& pts, vector<Element3D>& mesh, vector<array<double, 3>>& pts_b, vector<Element2D>& mesh_b, vector<int>& pid_loc)//need vtk file with point label
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3] >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7];
			for (int j = 0; j < 8; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]][0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]][1];
				mesh[i].pts[j][2] = pts[mesh[i].IEN[j]][2];
			}

		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	vector<array<int, 4>> face;
	vector<array<int, 4>> fc_sort;
	vector<array<int, 6>> mesh_fc(mesh.size());
	int fcloc[6][4] = { { 0, 3, 2, 1 },{ 0, 1, 5, 4 },{ 1, 2, 6, 5 },{ 2, 3, 7, 6 },{ 0, 4, 7, 3 },{ 4, 5, 6, 7 } };
	for (unsigned int i = 0; i < mesh.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			array<int, 4> tmp1 = { mesh[i].IEN[fcloc[j][0]], mesh[i].IEN[fcloc[j][1]], mesh[i].IEN[fcloc[j][2]], mesh[i].IEN[fcloc[j][3]] };
			array<int, 4> tmp2(tmp1);
			sort(tmp2.begin(), tmp2.end());
			vector<array<int, 4>>::iterator it = find(fc_sort.begin(), fc_sort.end(), tmp2);
			mesh_fc[i][j] = it - fc_sort.begin();
			if (it == fc_sort.end())
			{
				face.push_back(tmp1);
				fc_sort.push_back(tmp2);
			}
		}
	}
	vector<vector<int>> fc2mesh(face.size());
	for (unsigned int i = 0; i < mesh.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			fc2mesh[mesh_fc[i][j]].push_back(i);
		}
	}
	vector<int> pt_flag(pts.size(), 0);
	for (unsigned int i = 0; i < face.size(); i++)
	{
		if (fc2mesh[i].size() == 1)
		{
			Element2D etmp;
			for (int j = 0; j < 4; j++)
			{
				pt_flag[face[i][j]] = 1;
				etmp.IEN[j] = face[i][j];
				etmp.pts[j][0] = pts[face[i][j]][0];
				etmp.pts[j][1] = pts[face[i][j]][1];
				etmp.pts[j][2] = pts[face[i][j]][2];
			}
			mesh_b.push_back(etmp);
		}
	}
	int count(0);
	vector<int> pt_loc(pts.size(), -1);
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		if (pt_flag[i] == 1)
		{
			pts_b.push_back(pts[i]);
			pt_loc[i] = count++;
		}
	}
	for (unsigned int i = 0; i < mesh_b.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (pt_loc[mesh_b[i].IEN[j]] == -1)
			{
				cerr << "wrong index!\n"; getchar();
			}
			mesh_b[i].cnct[j] = pt_loc[mesh_b[i].IEN[j]];
		}
	}
	pid_loc = pt_loc;

	//cout << "# pts_b: " << pts_b.size() << "\n";
	//cout << "# mesh_b: " << mesh_b.size() << "\n";
	//getchar();
}

void ReadMeshLabel(string fn, vector<array<double, 3>>& pts, vector<int>& label, vector<Element3D>& mesh, vector<array<double, 3>>& pts_b, vector<Element2D>& mesh_b, vector<int>& pid_loc)//need vtk file with point label
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		label.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2] >>label[i]; 
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3] >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7];
			for (int j = 0; j < 8; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]][0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]][1];
				mesh[i].pts[j][2] = pts[mesh[i].IEN[j]][2];
			}

		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	vector<array<int, 4>> face;
	vector<array<int, 4>> fc_sort;
	vector<array<int, 6>> mesh_fc(mesh.size());
	int fcloc[6][4] = { { 0, 3, 2, 1 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 0, 4, 7, 3 }, { 4, 5, 6, 7 } };
	for (unsigned int i = 0; i < mesh.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			array<int, 4> tmp1 = { mesh[i].IEN[fcloc[j][0]], mesh[i].IEN[fcloc[j][1]], mesh[i].IEN[fcloc[j][2]], mesh[i].IEN[fcloc[j][3]] };
			array<int, 4> tmp2(tmp1);
			sort(tmp2.begin(), tmp2.end());
			vector<array<int, 4>>::iterator it = find(fc_sort.begin(), fc_sort.end(), tmp2);
			mesh_fc[i][j] = it - fc_sort.begin();
			if (it == fc_sort.end())
			{
				face.push_back(tmp1);
				fc_sort.push_back(tmp2);
			}
		}
	}
	vector<vector<int>> fc2mesh(face.size());
	for (unsigned int i = 0; i < mesh.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			fc2mesh[mesh_fc[i][j]].push_back(i);
		}
	}
	vector<int> pt_flag(pts.size(), 0);
	for (unsigned int i = 0; i < face.size(); i++)
	{
		if (fc2mesh[i].size() == 1)
		{
			Element2D etmp;
			for (int j = 0; j < 4; j++)
			{
				pt_flag[face[i][j]] = 1;
				etmp.IEN[j] = face[i][j];
				etmp.pts[j][0] = pts[face[i][j]][0];
				etmp.pts[j][1] = pts[face[i][j]][1];
				etmp.pts[j][2] = pts[face[i][j]][2];
			}
			mesh_b.push_back(etmp);
		}
	}
	int count(0);
	vector<int> pt_loc(pts.size(),-1);
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		if (pt_flag[i] == 1)
		{
			pts_b.push_back(pts[i]);
			pt_loc[i] = count++;
		}
	}
	for (unsigned int i = 0; i < mesh_b.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (pt_loc[mesh_b[i].IEN[j]] == -1)
			{
				cerr << "wrong index!\n"; getchar();
			}
			mesh_b[i].cnct[j] = pt_loc[mesh_b[i].IEN[j]];
		}
	}
	pid_loc = pt_loc;

	//cout << "# pts_b: " << pts_b.size() << "\n";
	//cout << "# mesh_b: " << mesh_b.size() << "\n";
	//getchar();
}

void ReadBezierElement(string fn, vector<Element3D>& mesh)
{

	string stmp;
	int npts, neles, nfunctions, itmp;

	string fname_cmat = fn + "_cmat.txt";
	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open())
	{
		fin_cmat >> neles;
		mesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin_cmat >> stmp >> nfunctions >> mesh[i].type;
			mesh[i].cmat.resize(nfunctions);
			mesh[i].IEN.resize(nfunctions);
			for (int j = 0; j < nfunctions; j++)
			{
				fin_cmat >> mesh[i].IEN[j];
			}
			for (int j = 0; j < nfunctions; j++)
			{
				for (int k = 0; k < 64; k++)
				{
					fin_cmat >> mesh[i].cmat[j][k];
				}

			}
		}
		fin_cmat.close();
		cout << "Bezier Matrices Loaded!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname_cmat << "!\n";
	}

	string fname_bzpt = fn + "_bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	if (fin_bzpt.is_open())
	{
		fin_bzpt >> npts;
		for (int e = 0; e < mesh.size(); e++)
		{
			mesh[e].pts.resize(bzpt_num);
			for (int i = 0; i < bzpt_num; i++)
			{
				fin_bzpt >> mesh[e].pts[i][0] >> mesh[e].pts[i][1] >> mesh[e].pts[i][2];
			}
		}
		fin_bzpt.close();
		cout << "Bezier Points Loaded!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname_bzpt << "!\n";
	}
}

void ReadVelocityField(string fn, vector<Element3D>& mesh)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < mesh.size(); i++)
		{
			fin >> mesh[i].velocity[0] >> mesh[i].velocity[1] >> mesh[i].velocity[2];
		}
		fin.close();
		cout << "Velocity Field Loaded!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void ReadVelocityFieldNode(string fn, vector<array<double, 3>>& pts, vector<array<double,3>>& velocity)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	velocity.resize(pts.size());
	if (fin.is_open())
	{
		for (int i = 0; i < pts.size(); i++)
		{
			fin >> velocity[i][0] >> velocity[i][1] >> velocity[i][2];
		}
		fin.close();
		cout << "Velocity Field Loaded!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void DataTrans2React(int ndof_b, const vector<int>& pid_loc, const vector<double>& CA, vector<double>& CA_b)
{
	CA_b.clear();
	CA_b.resize(ndof_b);
	for (unsigned int i = 0; i < pid_loc.size(); i++)
	{
		if (pid_loc[i] != -1)
		{
			CA_b[pid_loc[i]] = CA[i];
		}
	}
}

void DataTrans2Diffuse(int ndof, const vector<int>& pid_loc, const vector<double>& CA_b, vector<double>& CA)
{
	for (unsigned int i = 0; i < pid_loc.size(); i++)
	{
		if (pid_loc[i] != -1)
		{
			CA[i]=CA_b[pid_loc[i]] ;
		}
	}
}

void DataTrans2DiffuseBv(int ndof, const vector<int>& pid_loc, const vector<double>& Bv, vector<double>& Bv_all)
{
	Bv_all.clear();
	Bv_all.resize(ndof, 0.);
	for (unsigned int i = 0; i < pid_loc.size(); i++)
	{
		if (pid_loc[i] != -1)
		{
			Bv_all[i] = Bv[pid_loc[i]];
		}
	}
}