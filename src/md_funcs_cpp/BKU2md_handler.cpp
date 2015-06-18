#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>

#include "../utils_cpp/lib_includes.h"
#include "../classes/class_header.h"
#include "lib_md_includes.h"
#include "../math/jmath.h"
//#include "../utils_cpp/lib_various_tools.h"


using namespace std;

// ********************************************************************* //
// ******************memhandler CLASS MEMBER FUNCTIONS****************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Allocate MD Working Memory  ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::md_mem_alloc (MemHandler *data_mem,dataOutput* optfile)
{
	N = data_mem->N;	
	K = data_mem->K;	
	M = data_mem->Mtot;
	dt = data_mem->ipt_parms.dt;
	steps = data_mem->ipt_parms.steps;

	//Copy stating atom data
	for (int i = 0 ; i < N; ++i)
		{Adat.push_back(data_mem->atom_data[i].AtomMass());}

	//Produce bonding data stuff	
        for (int i = 0 ; i < K; ++i)
        {
		Bdat.push_back(data_mem->k[i]);
		Bdat[i].r0 = data_mem->r0[i];
	}

	optfile->ofile << "N: " << N << " K: " << K << " M: " << M << " dt: " << dt << "\n";
};

/*____________________________________________________________________________
                   ----Move Starting Vectors       ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::mv_starting_vectors (MemHandler *data_mem,dataOutput* optfile)
{
	//Copy initial position  and velocity vector from Memhandler 
        for (int i = 0 ; i < N; ++i)
        {
        	for (int j = 0 ; j < 3; ++j)
                {
			Adat[i].Pt.save(j,data_mem->pos_vec[j + i * 3]);
			Adat[i].Pi.save(j,data_mem->pos_vec[j + i * 3]);
			Adat[i].Vt.save(j,data_mem->vlc_vec[j + i * 3]);
			
		}
        }
        //Copy bonding data from Memhandler
        for (int i = 0; i < K; ++i)
        {
                Bdat[i].atom1 = data_mem->bonddata[i];
                Bdat[i].atom2 = data_mem->bonddata[i + K];
        }

	calc_bonds_per_atom(optfile);

        //Check if velocities need to be initialized
        double tMag = 0.0;
        for (int i = 0; i < N; ++i)
                {tMag += jsm::magnitude(Adat[i].Vt);}

        if (tMag < 1.0E-6)
        {
                optfile->ofile << "No starting velocities detected, running velocity initialization.\n";
                Velocity_Initialization(data_mem,optfile);
        } else {
                optfile->ofile << "Starting velocities detected, no need for velocity initialization.\n";
        }

        //Copy initial position  and velocity vector from Memhandler 
        for (int i = 0 ; i < N; ++i)
        {
                //Initialize the positions
                jsm::vec3<double> tm1 = Adat[i].Pt - Adat[i].Vt;
                jsm::vec3<double> tp1 = Adat[i].Pt + Adat[i].Vt;
                Adat[i].Ptm1 = tm1;
                Adat[i].Ptp1 = tp1;
        }

	//Print the position vectors
        optfile->ofile << "\nPositions tm1: \n";
        for (int i = 0; i < N; ++i)
                {optfile->ofile << " " << i << Adat[i].Ptm1 << "\n";}

	optfile->ofile << "\nPositions t: \n";
	for (int i = 0; i < N; ++i) 
		{optfile->ofile << " " << i << Adat[i].Pt << "\n";}

        optfile->ofile << "\nPositions tp1: \n";
        for (int i = 0; i < N; ++i)
                {optfile->ofile << " " << i << Adat[i].Ptp1 << "\n";}

        //Print the velocity vectors
        optfile->ofile << "\nInitial Velocities: \n";
        for (int i = 0; i < N; ++i)
                {optfile->ofile << " " << i << Adat[i].Vt << "\n";}

	optfile->ofile << "\nBonding: \n";
	for (int i = 0; i < K; ++i)
        { 
                optfile->ofile << " " << i << ": ";
                optfile->ofile << Bdat[i].atom1 << " <-> ";
                optfile->ofile << Bdat[i].atom2 << "\n";
        }

};

/*____________________________________________________________________________
                   ----Initialize Velocities       ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::Velocity_Initialization (MemHandler *data_mem,dataOutput* optfile)
{
        double kb = 1.3806488E-23;
	double Tr = data_mem->ipt_parms.temp;

	double Vscale = sqrt((3.0f * N * kb)/(double)M) * 1.0E+10 * dt;
	optfile->ofile << "\nVelocity Scale: " << Vscale << "\n";

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(-1.0,1.0);
	
	for (int i = 0; i < N; ++i)
	{
		optfile->ofile << "\n|----Calc velocities for atom: " << i << "----|\n";
		jsm::vec3<double> vf;
		for (int j = 0; j < Adat[i].NBonds; ++j)
		{
			int Abond = i;
			int Bbond = Adat[i].BA[j];

			jsm::vec3<double> va = Adat[Abond].Pt;
			jsm::vec3<double> vb = Adat[Bbond].Pt;
			jsm::vec3<double> vc = va - vb;
			jsm::vec3<double> vc2 = jsm::normalize(vc);

			optfile->ofile << "AtomA: " << Abond << " AtomB: " << Bbond << " Vc: " << vc2 << " Va: " << va << " Vb: " << vb << endl;

			vf += vc2;
			jsm::vec3<double> vfnorm = jsm::normalize(vf);
			vf = vfnorm;
		}

		
		double rN = distribution(generator);
		while (abs(rN) < 0.5 )
			{rN = distribution(generator);}

		jsm::vec3<double> vtemp = UniformScalar(vf,rN * Vscale);
		optfile->ofile << "vT: " << vtemp << endl;
		Adat[i].Vt = vtemp;

		optfile->ofile << "|--------------------------------|\n\n";
	}

	optfile->ofile << "\n|-------Scaling Bond Momentum-------|\n";
	for (int i = 0; i < K; ++i)
	{
		jsm::vec3<double> aP; //Calculate average linear momentum of the bond
			
		int atom1 = Bdat[i].atom1;
		int atom2 = Bdat[i].atom2;

		jsm::vec3<double> wv1 = Adat[atom1].Vt;
		jsm::vec3<double> wv2 = Adat[atom2].Vt;

		jsm::vec3<double> swv1 = jsm::UniformScalar(wv1,Adat[atom1].AM);
		jsm::vec3<double> swv2 = jsm::UniformScalar(wv2,Adat[atom2].AM);

		aP += swv1;
		aP += swv2;

		jsm::vec3<double> aPa = jsm::UniformScalar(aP,0.5);
		optfile->ofile << "Bond(" << i << "): Total Bond Momentum Vector (Before Correction): " << aPa << "\n";
	
		jsm::vec3<double> iP1 = jsm::UniformScalar(Adat[atom1].Vt,Adat[atom1].AM);
		jsm::vec3<double> iP2 = jsm::UniformScalar(Adat[atom2].Vt,Adat[atom2].AM);

		jsm::vec3<double> iPi1 = iP1 - aPa;
		jsm::vec3<double> iPi2 = iP2 - aPa;

		jsm::vec3<double> iVi1 = UniformScalar(iPi1, 1/(double)Adat[atom1].AM);
		jsm::vec3<double> iVi2 = UniformScalar(iPi2, 1/(double)Adat[atom2].AM);

		Adat[atom1].Vt = iVi1;
		Adat[atom2].Vt = iVi2;

		optfile->ofile << "\n|------------------------------|\n";
	}

	for (int i = 0; i < N; ++i)
        {
		double delta = 1.0E-10;
		if (abs(Adat[i].Vt.x) < delta)
			{Adat[i].Vt.x = 0.0f;}
                if (abs(Adat[i].Vt.y) < delta)
                        {Adat[i].Vt.y = 0.0f;}
                if (abs(Adat[i].Vt.z) < delta)
                        {Adat[i].Vt.z = 0.0f;}
	}
};

/*____________________________________________________________________________
                   ----Input CM Shift              ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::CM_Calc (MemHandler *data_mem,dataOutput* optfile)
{
        double num[3],CM[3];
        num[0] = 0.00; num[1] = 0.00; num[2] = 0.00;

        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < 3; ++j)
                	{num[j] += Adat[i].AM * Adat[i].Pt.fetch(j);}
        }

        for (int j = 0; j < 3; ++j)
        	{CM[j] = num[j] / (double) M;}

        //optfile->ofile << "N: " << N << "\n";
        optfile->ofile << "CM[x]: " << CM[0] << " CM[y]: " << CM[1] <<" CM[z]: " << CM[2] <<"\n";
}

/*____________________________________________________________________________
                   ----Determine Bonds Per Atom    ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::calc_bonds_per_atom(dataOutput* optfile)
{
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < K; ++j)
                {
                        if (i == Bdat[j].atom1)
                        {
				Adat[i].BA.push_back(Bdat[j].atom2);
				Adat[i].Kc.push_back(Bdat[j].Kc);
				Adat[i].r0.push_back(Bdat[j].r0);
				++Adat[i].NBonds;
			}

                        if (i == Bdat[j].atom2)
                        {
				Adat[i].BA.push_back(Bdat[j].atom1);
                                Adat[i].Kc.push_back(Bdat[j].Kc);
                                Adat[i].r0.push_back(Bdat[j].r0);
				++Adat[i].NBonds;
			}
                }
        }

        //for (int i = 0; i < N; ++i)
        //        {optfile->ofile << "(" << i << "): " << bonds[i] << "\n";}
}

/*____________________________________________________________________________
                   ----Calculate Forces            ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::calc_forces(MemHandler *data_mem,dataOutput* optfile)
{
	//Calculate force vector for each atom.
	#pragma omp parallel for shared(data_mem)
        for (int i = 0; i < N; ++i)
        {
		//Determine Total Force Components (for each F[x] = F[0],F[y] = F[1],F[z] = F[2])
		jsm::vec3<double> Ftot;

		//Calculate forces for each atom
		for (int l = 0; l < Adat[i].NBonds; ++l)
		{
			int atomA = i; // Atom A
			int atomB = Adat[i].BA[l]; //Atom B

			jsm::vec3<double> pA = Adat[atomA].Pt; 
			jsm::vec3<double> pB = Adat[atomB].Pt;

			double Kc = Adat[atomA].Kc[l];
			double r0 = Adat[atomA].r0[l];
		
                        jsm::vec3<double> x = pA - pB;
                        jsm::vec3<double> x0 = jsm::normalize(x) * r0;

			jsm::vec3<double> rv = x - x0;

			cout << "TEST\n";

			if (abs(rv.x) < 1.0E-13)
				{rv.x = 0.0;}

                        if (abs(rv.y) < 1.0E-13)double Tr = data_mem->ipt_parms.temp;
                                {rv.y = 0.0;}

                        if (abs(rv.z) < 1.0E-13)
                                {rv.z = 0.0;}

			optfile->ofile << "Kc: " << Kc << " x: " << x << " x0: " << x0 << " x-x0: " << rv << "\n";

			jsm::vec3<double> Fvec = (rv) * (-Kc);
			Ftot += Fvec;
		}

		Adat[i].Ft = Ftot;
	}

	optfile->ofile << "Force Vectors: \n";
	for (int l = 0; l < N; ++l)
		{optfile->ofile << "Force (" << l << ")" << Adat[l].Ft << "\n";}
	optfile->ofile << "\n";
}

/*____________________________________________________________________________
                   ----Verlet Integration          ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::verlet_integration(MemHandler *data_mem,dataOutput* optfile)
{
	double dt2 = (double)(dt * dt);

        for (int i = 0; i < N; ++i)
	{
		double AM = (double)Adat[i].AM;

		jsm::vec3<double> av = Adat[i].Ft * (dt2/(double)AM);
		jsm::vec3<double> t1 = Adat[i].Pt * 2.0;

                if (abs(av.x) < 1.0E-13)
                        {av.x = 0.0;}

                if (abs(av.y) < 1.0E-13)
                        {av.y = 0.0;}

                if (abs(av.z) < 1.0E-13)
                        {av.z = 0.0;}

		optfile->ofile << "T1: " << t1 << " T2: " << Adat[i].Ptm1 << " T3: " << av << "\n";
		Adat[i].Ptp1 = t1 - Adat[i].Ptm1 + av;
	}
}

/*____________________________________________________________________________
                   ----Shift Pos Vectors           ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 12/16/2014   ----
                   ----Modified By:                ----
*/
void md_funcs::shift_pos_vec(dataOutput* optfile)
{
	for (int i = 0; i < N; ++i)
        {
                optfile->ofile << "ATOM 1 POSITIONS: \n";
                optfile->ofile << "Pos(t-1) - x: " << Adat[i].Ptm1 << "\n";
                optfile->ofile << "Pos(t) - x: " << Adat[i].Pt << "\n";
                optfile->ofile << "Pos(t+1) - x: " << Adat[i].Ptp1 << "\n";
        

		Adat[i].Ptm1 = Adat[i].Pt;
		Adat[i].Pt = Adat[i].Ptp1;
	}
}

/*____________________________________________________________________________
                   ----Calculate Potential       ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
double md_funcs::calc_potential(int bond,dataOutput* optfile)
{
	double r_T = jsm::magnitude(Adat[Bdat[bond].atom1].Pt - Adat[Bdat[bond].atom2].Pt);
        float Kc = Bdat[bond].Kc;
        float r0 = Bdat[bond].r0;

	return Kc * (r_T - r0) * (r_T - r0);
}

/*____________________________________________________________________________
                   ----Calculate Kinetic         ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
double md_funcs::calc_kinetic(int atom,dataOutput* optfile)
{
        double Tdist = jsm::magnitude(Adat[atom].Ptm1 - Adat[atom].Ptp1);
	double AM = (double)Adat[atom].AM;

	double v = Tdist / (double)(2.0 * dt);

	optfile->ofile << "Atom: " << atom << " Velocity: " << v << "\n";

        return AM * v * v;
}

/*____________________________________________________________________________
                   ----Calculate E Total         ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
double md_funcs::calc_E_total(MemHandler *data_mem,dataOutput* optfile)
{
	//Calculate the total potential energy stored in each bond
	double Vtot = 0;
	for (int i = 0; i < K; ++i)
		{Vtot += calc_potential(i,optfile);}
	Vtot = 0.5 * Vtot;

	//Calculate the total kinetic energy of the atomic motions
        double Ktot = 0;
        for (int i = 0; i < N; ++i)
        	{Ktot += calc_kinetic(i,optfile);}
        Ktot = 0.5 * Ktot;

	//Ev_tp1 = Ev_t;
	//Ev_t = Vtot;

	double kb = 0.0019872041;
	double T = (2 * Ktot) / (double)(3 * N * kb);

	//SCALE VELOCITIES
	scale_velocities(data_mem,optfile,T);
	Tavg = ((Tavg * (step - 1)) + T)/(double)(step);
	optfile->ofile << "T(avg): " << Tavg << " T(inst): " << T << "\n";

	double Etot = Ktot + Vtot;

	optfile->ofile << "Etot: " << Etot << " Vtot: " << Vtot << " Ktot: " << Ktot << "\n";

        ofstream graph[4];
	stringstream os[4];

	os[0] << data_mem->ipt_parms.data_dir << "Etot_graph.dat";
	os[1] << data_mem->ipt_parms.data_dir << "Vtot_graph.dat";
	os[2] << data_mem->ipt_parms.data_dir << "Ktot_graph.dat";
	os[3] << data_mem->ipt_parms.data_dir << "Ttot_graph.dat";
        graph[0].open(os[0].str().c_str(),ios::out | ios::app);
        graph[1].open(os[1].str().c_str(),ios::out | ios::app);
        graph[2].open(os[2].str().c_str(),ios::out | ios::app);
        graph[3].open(os[3].str().c_str(),ios::out | ios::app);

        int atom1 = Bdat[0].atom1;
        int atom2 = Bdat[0].atom2;

	float radius = jsm::magnitude(Adat[atom1].Pt - Adat[atom2].Pt); 

        graph[0] << radius << "  " << Etot << "\n";
        graph[1] << radius << "  " << Vtot << "\n";
        graph[2] << radius << "  " << Ktot << "\n";
        graph[3] << step * dt << "  " << T << "\n";

        graph[0].close();
        graph[1].close();
        graph[2].close();
        graph[3].close();

	return Etot;	
};

/*____________________________________________________________________________
                   ----Produce MD Output         ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
void md_funcs::produce_md_out(MemHandler *data_mem,dataOutput* optfile)
{
        ofstream mdout;
        stringstream os;
	os << data_mem->ipt_parms.data_dir << "mdout.xyz";
        mdout.open(os.str().c_str(),ios::out | ios::app);

	mdout << N <<"\n\n";

	for (int i = 0; i < N; ++i)
	{
 	       mdout << data_mem->atom_data[i].AtomLetter() << "    " << Adat[i].Pt.x << "      " << Adat[i].Pt.y << "      " << Adat[i].Pt.z << "\n";
	}

        mdout.close();
};

/*____________________________________________________________________________
                   ----Produce Restart File      ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 1/29/2014 ----
                   ----Modified By:              ----
*/
/*void md_funcs::produce_restart(MemHandler *data_mem,dataOutput* optfile)
{
        ofstream mdout;
        stringstream os;
        os << data_mem->ipt_parms.data_dir << "restart.rst";
        mdout.open(os.str().c_str(),ios::out | ios::app);

        mdout << N <<"\n\n";

        for (int i = 0; i < N; ++i)
        {
               mdout << data_mem->atom_data[i].AtomLetter() << "    " << Adat[i].Pt.x << "      " << Adat[i].Pt.y << "      " << Adat[i].Pt.z << "\n";
        }

        mdout.close();
};*/

/*____________________________________________________________________________
                   ----Calculate RMSD            ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 1/21/2014  ----
                   ----Modified By:              ----
*/
void md_funcs::calculate_rmsd(MemHandler *data_mem,dataOutput* optfile)
{
	++rmsdCount;
	float rmsd = 0;
	for (int i = 0; i < N; ++i)
	{
		float d = jsm::magnitude(Adat[i].Pt - Adat[i].Pi);
		rmsd += d * d;
	}

	rmsd = sqrt(rmsd / (float)N);
	rmsdSum += rmsd;

	//Save RSMD Graphs
        ofstream graph[2];
        stringstream os[2];

        os[0] << data_mem->ipt_parms.data_dir << "AvgRMSD.dat";
        os[1] << data_mem->ipt_parms.data_dir << "InstRMSD.dat";

        graph[0].open(os[0].str().c_str(),ios::out | ios::app);
        graph[1].open(os[1].str().c_str(),ios::out | ios::app);

        graph[0] << step * dt << "  " << rmsdSum / (float)rmsdCount << "\n";
        graph[1] << step * dt << "  " << rmsd << "\n";

        graph[0].close();
        graph[1].close();
};

/*____________________________________________________________________________
                   ----Velocity Scaling          ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 1/21/2014  ----
                   ----Modified By:              ----
*/
void md_funcs::scale_velocities(MemHandler *data_mem,dataOutput* optfile,double T)
{
	int HS = data_mem->ipt_parms.heatSteps;
        double Tr = data_mem->ipt_parms.temp;
	double tauM = data_mem->ipt_parms.tau;

        //std::default_random_engine generator;
        //std::uniform_real_distribution<double> distribution(-0.01,0.01);

	if (step < HS)
	{
		Tr = ceil(Tr * (step / (double)HS));
	}

	double tau = tauM * dt;

	double lambda = sqrt(1 + (dt/(double)tau)*((Tr/(double)T) - 1));

	if (step % 100 == 0)
	{
		cout << "Step: " << step << " LAMBDA: " << lambda << " dt/tau: " << dt/(double)tau << " SetTemp: " << Tr << " Temp: " << T << "\n";
	}

	for (int i = 0; i < N; ++i)
        {
		//double rN = 1;
		//double rN = distribution(generator);

		int atom1 = i;
        	jsm::vec3<double> vel_s1 = jsm::UniformScalar(Adat[atom1].Ptp1 - Adat[atom1].Pt,lambda);
        	Adat[atom1].Ptp1 = Adat[atom1].Pt + vel_s1;
	}
};


