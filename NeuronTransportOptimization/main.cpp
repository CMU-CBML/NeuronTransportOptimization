#include <iostream>
#include <vector>
#include <array>
#include "Laplace.h"
#include "diff_react.h"
#include "Test.h"
#include <sstream>
#include <iomanip>
using namespace std;

int main()
{
	double tstep(0.5);
	int nstep(11);
	double CAi, CAs,length;//for test
	double error;//for test
	double time(0.);//for test
	vector<double> var;
	vector <double> Error;
	vector<array<double, 3>> cpts;
	vector<array<double, 3>> cpts_b;
	vector<array<double, 3>> velocity_node;
	vector<Element3D> tmesh;
	vector<Element2D> tmesh_b;
	vector<Element3D> bzmesh;
	vector<int>label;
	vector<int> pid_loc;
	vector<double> CA0, N_plus0, N_minus0, Bv, Bv_all;
	vector<double> CA_temp, N_plus_temp, N_minus_temp;
	SetVariables(var);

	////Cylinder Model
	//ReadMeshLabel("../io/IGA/cylinder/cylinder_controlmesh.vtk", cpts, label, tmesh, cpts_b, tmesh_b, pid_loc);
	//ReadBezierElement("../io/IGA/cylinder/cylinder", bzmesh);
	//ReadVelocityFieldNode("../io/IGA/cylinder/hexmesh_node_velocityfield.txt", cpts, velocity_node);
	//string fn("../io/IGA/cylinder/CA1.0_NP2.0_v1.0_k0.5_dt0.1/cylinder_hex_");

	//Cylinder Photoactivation Model
	ReadMeshLabel("../io/IGA/cylinder_photoactivation/cylinder_controlmesh.vtk", cpts, label, tmesh, cpts_b, tmesh_b, pid_loc);
	ReadBezierElement("../io/IGA/cylinder_photoactivation/cylinder", bzmesh);
	ReadVelocityFieldNode("../io/IGA/cylinder_photoactivation/hexmesh_node_velocityfield.txt", cpts, velocity_node);
	string fn("../io/IGA/cylinder_photoactivation/CA1.0_NP2.0_v0.5_k0.5_dt0.5/cylinder_hex_");

	////Bifurcation Smooth Model
	//ReadMeshLabel("../io/IGA/bifurcation_smooth1/bifurcation_smooth_controlmesh.vtk", cpts, label, tmesh, cpts_b, tmesh_b, pid_loc);
	//ReadBezierElement("../io/IGA/bifurcation_smooth1/bifurcation_smooth", bzmesh);
	//ReadVelocityFieldNode("../io/IGA/bifurcation_smooth1/hexmesh_node_velocityfield.txt", cpts, velocity_node);
	//string fn("../io/IGA/bifurcation_smooth1/CA1.0_NP2.0_v1.0_k0.1_dt0.1/bifurcation_smooth_");

	////Bifurcation Photoactivation Model
	//ReadMeshLabel("../io/IGA/bifurcation_photoactivation/bifurcation_smooth_controlmesh.vtk", cpts, label, tmesh, cpts_b, tmesh_b, pid_loc);
	//ReadBezierElement("../io/IGA/bifurcation_photoactivation/bifurcation_smooth", bzmesh);
	//ReadVelocityFieldNode("../io/IGA/bifurcation_photoactivation/hexmesh_node_velocityfield.txt", cpts, velocity_node);
	//string fn("../io/IGA/bifurcation_photoactivation/CA1.0_NP2.0_v1.0_k0.5_dt1.0/bifurcation_smooth_");

	////3 Bifurcations Neuron Model
	//ReadMeshLabel("../io/IGA/3bifurcation/3bifurcation_controlmesh.vtk", cpts, label, tmesh, cpts_b, tmesh_b, pid_loc);
	//ReadBezierElement("../io/IGA/3bifurcation/3bifurcation", bzmesh);
	//ReadVelocityFieldNode("../io/IGA/3bifurcation/hexmesh_node_velocityfield.txt", cpts, velocity_node);
	//string fn("../io/IGA/3bifurcation/CA1.0_NP2.0_v0.2_k0.1_dt0.01/3bifurcation_hex_");


	////3 Bifurcations Neuron Smooth Model
	//ReadMeshLabel("../io/neuron_3bifurcations_smooth/neuron_label_smooth.vtk", pts, label, mesh, pts_b, mesh_b, pid_loc);
	//ReadVelocityFieldNode("../io/neuron_3bifurcations_smooth/hexmesh_node_velocityfield.txt", pts, velocity_node);
	//string fn("../io/neuron_3bifurcations_smooth/CA1.0_NP2.0_v1.0_k0.1_dt0.01/neuron_hex_");

	SetInitialCondition(cpts.size(), cpts_b.size(), CA0, N_plus0, N_minus0, cpts, label, pid_loc);
	SetTempData(cpts.size(), cpts_b.size(), CA_temp, N_plus_temp, N_minus_temp);

	Laplace diffuse;
	diffuse.InitializeProblem(bzmesh, velocity_node, CA0, N_plus0, N_minus0, var, tstep); //CA = CA0
	//diffuse.Run(pts, velocity_node, mesh, label, pid_loc);

	//TEST test1;
	//test1.InitializeProblem(CA0, var, CAi, CAs, tstep, length);

	for (int i = 0; i < nstep; ++i)
	{
		stringstream ss;
		ss << i;

		cout << "Step: " << i << "\n";
		time += tstep;
		CA_temp = diffuse.GetCA();
		N_plus_temp = diffuse.GetNplus();
		N_minus_temp = diffuse.GetNminus();
		
		//if (i==0 ||(i+1) % 10 == 0)
		{
			//diffuse.VisualizeVTK(cpts, tmesh, fn + ss.str());
			//diffuse.VisualizeVTK_1(bzmesh, fn + ss.str());
			diffuse.VisualizeVTK_2(bzmesh, fn + ss.str());
			//test1.VisualizeVTK(pts, mesh, fn + ss.str());
			//test1.VisualizeVTK_ERROR(pts, mesh, fn + ss.str());
		}

		diffuse.Run(cpts, velocity_node, tmesh, bzmesh, label, pid_loc);
		cout << "N0 Error" << MaxDifference(CA_temp, diffuse.GetCA()) << "\n";
		cout << "Nplus Error" << MaxDifference(N_plus_temp, diffuse.GetNplus()) << "\n";

		//if (i==0 ||(i+1) % 10 == 0)
		{
			//diffuse.VisualizeVTK(cpts, tmesh, fn + ss.str());
			//diffuse.VisualizeVTK_1(bzmesh, fn + ss.str());
			//diffuse.VisualizeVTK_2(bzmesh, fn + ss.str());
			//test1.VisualizeVTK(pts, mesh, fn + ss.str());
			//test1.VisualizeVTK_ERROR(pts, mesh, fn + ss.str());
		}


		//string fn44;
		//fn44 = "../io/test6/L2Error_Adaptive.txt";
		//ofstream fout44;
		//fout44.open(fn44, ios::app);
		//fout44 << time << " " << error<< "\n";
		////fout44 << time << " " << error << " " << MaxDifference(CA_temp, diffuse.CA) << " " << MaxDifference(test1.CA , diffuse.CA) << "\n";
		//cout << "FinalError_L2 " << error << "\n";
		//cout << "MaxError" << MaxDifference(diffuse.CA, test1.CA) << "\n";
		cout << "Step " << i << " done!\n";
		//if (MaxDifference(CA_temp,diffuse.CA)<1.e-7)//don't have analytical solution
		//{
		//	break;
		//}
	}	
	

	cout<<"DONE!\n";
	getchar();
	return 0;
}