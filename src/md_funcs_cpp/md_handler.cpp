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
    {
        Adat.push_back(data_mem->atom_data[i].AtomMass());
    }

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
    {
        tMag += jsm::magnitude(Adat[i].Vt);
    }

    if (tMag < 1.0E-6)
    //if (true)
    {
        optfile->ofile << "No starting velocities detected, running velocity initialization.\n";
        //cout << "No starting velocities detected, running velocity initialization.\n";
        Velocity_Initialization(data_mem,optfile);
    }
    else
    {
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
    optfile->ofile << "\nPositions: \n";
    for (int i = 0; i < N; ++i)
    {
        optfile->ofile << " tm1: " << i << Adat[i].Ptm1 << "\n";
        optfile->ofile << " t: " << i << Adat[i].Pt << "\n";
        optfile->ofile << " tp1: " << i << Adat[i].Ptp1 << "\n";
    }

    //Print the velocity vectors
    optfile->ofile << "\nInitial Velocities: \n";
    for (int i = 0; i < N; ++i)
    {
        optfile->ofile << " " << i << Adat[i].Vt << "\n";
    }

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
                   ----Date Modified: 3/6/2015   ----
                   ----Modified By:                ----
*/
void md_funcs::Velocity_Initialization (MemHandler *data_mem,dataOutput* optfile)
{
    double kb = 1.3806488E-23;
    //double Tr = data_mem->ipt_parms.temp;

    //double Vscale = sqrt((3.0f * N * kb)/(double)M) * 1.0E+10 * dt;
    double Vscale = sqrt((N * kb)/(double)M) * 1.0E+10 * dt;
    optfile->ofile << "\nVelocity Scale: " << Vscale << "\n";
    cout << "\nVelocity Scale: " << Vscale << "\n";

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.00001,-0.00001);

    for (int i = 0; i < N; ++i)
    {
        optfile->ofile << "\n|----Calc velocities for atom: " << i << "----|\n";
        //jsm::vec3<double> vf;
        /*for (int j = 0; j < Adat[i].NBonds; ++j)
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
        }*/


        //double rN = distribution(generator);
        //while (abs(rN) < 0.5 )
        //        {rN = distribution(generator);}

        //optfile->ofile << "rN: " << rN << endl;
        //jsm::vec3<double> vtemp = UniformScalar(vf,rN * Vscale);
        //optfile->ofile << "vf: " << vf << endl;
        jsm::vec3<double> vtemp(0.0,0.0,pow(-1,i) * Vscale);
        //optfile->ofile << "vT: " << vtemp << endl;
        Adat[i].Vt = vtemp;
        optfile->ofile << "vT(" << i << "): " << Adat[i].Vt << endl;
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
    }

    for (int i = 0; i < N; ++i)
    {
        double delta = 1.0E-10;
        if (abs(Adat[i].Vt.x) < delta)
        {
            Adat[i].Vt.x = 0.0f;
        }
        if (abs(Adat[i].Vt.y) < delta)
        {
            Adat[i].Vt.y = 0.0f;
        }
        if (abs(Adat[i].Vt.z) < delta)
        {
            Adat[i].Vt.z = 0.0f;
        }
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
    num[0] = 0.00;
    num[1] = 0.00;
    num[2] = 0.00;

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            num[j] += Adat[i].AM * Adat[i].Pt.fetch(j);
        }
    }

    for (int j = 0; j < 3; ++j)
    {
        CM[j] = num[j] / (double) M;
    }



    //optfile->ofile << "N: " << N << "\n";
    optfile->ofile << "CM[x]: " << CM[0] << " CM[y]: " << CM[1] <<" CM[z]: " << CM[2] <<"\n";

    if (step % 20 == 0)
    {
        optfile->graph[7] << step * dt * 1.0E12 << "  " << CM[2] << "\n";
    }
};

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

    /*optfile->ofile << "|************SAVED BONDS************|" << endl;
        for (int i = 0; i < N; ++i)
        {
        	for (int j = 0; j < (int)Adat[i].BA.size(); ++j)
    	{
    		optfile->ofile << "Atom(" << i << ") BONDED: " << Adat[i].BA[j] << "\n";
    	}
    }*/
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
    //double stime=omp_get_wtime();

    //#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        //Determine Total Force Components (for each F[x] = F[0],F[y] = F[1],F[z] = F[2])
        jsm::vec3<double> Ftot;

        int atomA = i; // Atom A
        jsm::vec3<double> pA = Adat[atomA].Pt;

        //Calculate forces for each atom
        for (int l = 0; l < Adat[i].NBonds; ++l)
        {
            int atomB = Adat[i].BA[l];
            double Kc = Adat[atomA].Kc[l];
            double r0 = Adat[atomA].r0[l];
            jsm::vec3<double> pB = Adat[atomB].Pt;

            jsm::vec3<double> x = pA - pB;
            jsm::vec3<double> x0 = jsm::normalize(x) * r0;

            jsm::vec3<double> rv = x - x0;

            //cout << "TEST\n";

            if (abs(rv.x) < 1.0E-14)
            {
                rv.x = 0.0;
            }

            if (abs(rv.y) < 1.0E-14)
            {
                rv.y = 0.0;
            }

            if (abs(rv.z) < 1.0E-14)
            {
                rv.z = 0.0;
            }

            //optfile->ofile << "Kc: " << Kc << " x: " << x << " x0: " << x0 << " x-x0: " << rv << "\n";

            jsm::vec3<double> Fvec = (rv) * (-Kc);
            Ftot += Fvec;
        }

        Adat[i].Ft = Ftot;
    }


    //double etime=omp_get_wtime();

    //cout << "Force Time: " << etime-stime << endl;
    //optfile->ofile << "Force Vectors: \n";
    //for (int l = 0; l < N; ++l)
    //	{optfile->ofile << "Force (" << l << ")" << Adat[l].Ft << "\n";}
    //optfile->ofile << "\n";
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
        //optfile->ofile << "Atom: " << i << " Ptm1: " << Adat[i].Ptm1 << " Pt: " << Adat[i].Pt << " Ptp1: " << Adat[i].Ptp1 << "\n";
        double AM = (double)Adat[i].AM;

        jsm::vec3<double> av = Adat[i].Ft * (dt2/(double)AM);
        jsm::vec3<double> t1 = Adat[i].Pt * 2.0;

        if (abs(av.x) < 1.0E-7)
        {
            av.x = 0.0;
        }

        if (abs(av.y) < 1.0E-7)
        {
            av.y = 0.0;
        }

        if (abs(av.z) < 1.0E-7)
        {
            av.z = 0.0;
        }

        double x = Adat[i].Pt.x;
        double y = Adat[i].Pt.y;
        double z = t1.z - Adat[i].Ptm1.z + av.z;
        jsm::vec3<double> Ptp1(x,y,z);
        Adat[i].Ptp1 = Ptp1;
        //optfile->ofile << "Atom: " << i << " Ptm1: " << Adat[i].Ptm1 << " Pt: " << Adat[i].Pt << " Ptp1: " << Adat[i].Ptp1 << " av: " << av << "\n";
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
        //optfile->ofile << "ATOM 1 POSITIONS: \n";
        //optfile->ofile << "Pos(t-1) - x: " << Adat[i].Ptm1 << "\n";
        //optfile->ofile << "Pos(t) - x: " << Adat[i].Pt << "\n";
        //optfile->ofile << "Pos(t+1) - x: " << Adat[i].Ptp1 << "\n";


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

    //optfile->ofile << "Atom: " << atom << " Velocity: " << v << " TDist: " << Tdist << " AM: " << AM <<  "\n";

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
    {
        Vtot += calc_potential(i,optfile);
    }
    Vtot = 0.5 * Vtot;

    //Calculate the total kinetic energy of the atomic motions
    double Ktot = 0;
    for (int i = 0; i < N; ++i)
    {
        Ktot += calc_kinetic(i,optfile);
    }
    Ktot = 0.5 * Ktot;

    //Ev_tp1 = Ev_t;
    //Ev_t = Vtot;

    double kb = 0.0019872041;
    //double T = (2 * Ktot) / (double)(3 * N * kb);
    double T = (2 * Ktot) / (double)(N * kb);

    //SCALE VELOCITIES
    scale_velocities(data_mem,optfile,T);
    Tavg = ((Tavg * (step - 1)) + T)/(double)(step);
    optfile->ofile << "T(avg): " << Tavg << " T(inst): " << T << "\n";

    double Etot = Ktot + Vtot;

    optfile->ofile << "Etot: " << Etot << " Vtot: " << Vtot << " Ktot: " << Ktot << "\n";

    if (step % 50 == 0)
    {
        int atom1 = Bdat[0].atom1;
        int atom2 = Bdat[0].atom2;

        float radius = jsm::magnitude(Adat[atom1].Pt - Adat[atom2].Pt);

        //cout << "SHOULD BE SAVING DATA!\n";

        optfile->graph[0] << radius << "  " << Etot << "\n";
        optfile->graph[1] << radius << "  " << Vtot << "\n";
        optfile->graph[2] << radius << "  " << Ktot << "\n";
        optfile->graph[3] << step * dt << "  " << T << "\n";
    }

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
    optfile->graph[6] << N <<"\n\n";

    for (int i = 0; i < N; ++i)
    {
        optfile->graph[6] << data_mem->atom_data[i].AtomLetter() << "    " << Adat[i].Pt.x << "      " << Adat[i].Pt.y << "      " << Adat[i].Pt.z << "\n";
    }
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

    vector<jsm::vec3<double>> position;
    position.resize(N);
    ShiftCM(position,optfile);

    ++rmsdCount;
    float rmsd = 0;
    for (int i = 0; i < N; ++i)
    {
        float d = position[i].z - Adat[i].Pi.z;
        rmsd += d * d;
    }

    rmsd = sqrt(rmsd / (float)N);
    rmsdSum += rmsd;

    //Save data after heating for statistics
    int HS = data_mem->ipt_parms.heatSteps;
    if (step > HS)
    {
        RMSD.push_back(rmsd);
    }

    //Save RSMD Graphs
    if(step % 50 == 0)
    {
        optfile->graph[4] << step * dt * 1.0E12 << "  " << rmsdSum / (float)rmsdCount << "\n";
        optfile->graph[5] << step * dt * 1.0E12 << "  " << rmsd << "\n";
    }
};

/*____________________________________________________________________________
                   ----Calculate RMSD            ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 1/21/2014  ----
                   ----Modified By:              ----
*/
void md_funcs::ShiftCM(vector<jsm::vec3<double>> &position,dataOutput* optfile)
{
    optfile->ofile << "\nCalculating RMSD CM Shift\n";
    double Ztotal = 0;
    for (int i = 0; i < N; ++i)
    {
        Ztotal += Adat[i].Pt.z * Adat[i].AM;
    }
    double Zcm = Ztotal / (double)M;

    optfile->ofile << "Zcm: " << Zcm << endl;

    for (int i = 0; i < N; ++i)
    {
        position[i].z = Adat[i].Pt.z - Zcm;
        //optfile->ofile << "Initial Pos: " << Adat[i].Pt << " Shifted Pos: " << position[i] << endl;

        if (step % 500 == 0)
        {
            Adat[i].Ptm1.z = Adat[i].Ptm1.z - Zcm;
            Adat[i].Pt.z = Adat[i].Pt.z - Zcm;
            Adat[i].Ptp1.z = Adat[i].Ptp1.z - Zcm;
        }
    }
    optfile->ofile << "\n";
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

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.00001,0.00001);

    std::default_random_engine generator2;
    std::uniform_real_distribution<double> distribution2(1.0,0.0);

    if (step < HS)
    {
        Tr = ceil(Tr * (step / (double)HS));
    }

    double tau = tauM * dt;
    //double tau = tauM;

    //cout << "Tr: " << Tr << endl;
    double lambda = sqrt(1 + (dt/(double)tau)*((Tr/(double)T) - 1));

    if (step % 100 == 0)
    {
        cout << "Step: " << step << " LAMBDA: " << lambda << " dt/tau: " << dt/(double)tau << " SetTemp: " << Tr << " Temp: " << T << "\n";
    }

    vector<jsm::vec3<double>> vel_scales;
    vel_scales.resize(N);

    if (!data_mem->coupled)
    {
        for (int i = 0; i < K; ++i)
        {
            if (distribution(generator2) < 0.5)
            {
                //double rN = 0;
                double rN = distribution(generator);
                //int atom1 = i;
                int atom1 = Bdat[i].atom1;
                int atom2 = Bdat[i].atom2;

                jsm::vec3<double> vel_s1 = jsm::zScalar(Adat[atom1].Ptp1 - Adat[atom1].Ptm1,lambda/2.0f);

                jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1.z+rN);
                //optfile->ofile << "SCALING VELOCITY: " << vel_sa1 << endl;
                //jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1.z);

                //Adat[atom1].Ptp1 = Adat[atom1].Pt + vel_sa1;
                //Adat[atom2].Ptp1 = Adat[atom2].Pt - vel_sa1;
                vel_scales[atom1] += vel_sa1;
                vel_scales[atom2] -= vel_sa1;

            }
        }
    }
    else
    {
        double sumRandV=0;
        vector<double> rN;

        for (int i=0; i<N; ++i)
        {
            rN.push_back(distribution(generator));
            sumRandV += rN[i];
        }

        sumRandV /= (double)N;
        for (int i=0; i<N; ++i)
        {
            rN[i] = rN[i] - sumRandV;
        }

        for (int i = 0; i < N; ++i)
        {
            //if (distribution(generator2) < 0.5)
            //{
            //double rN = 0;
            //double rN = distribution(generator);
            int atom1 = i;
            //int atom1 = Bdat[i].atom1;
            //int atom2 = Bdat[i].atom2;

            //optfile->ofile << "Atom: " << atom1 << " Ptm1: " << Adat[atom1].Ptm1 << " Pt: " << Adat[atom1].Pt << " Ptp1: " << Adat[atom1].Ptp1 <<"\n";
            //optfile->ofile << "Atom: " << atom2 << " Ptm1: " << Adat[atom2].Ptm1 << " Pt: " << Adat[atom2].Pt << " Ptp1: " << Adat[atom2].Ptp1 <<"\n";
            jsm::vec3<double> vel_s1 = jsm::zScalar(Adat[atom1].Ptp1 - Adat[atom1].Ptm1,lambda/2.0f);

            jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1.z+rN[i]);
            //optfile->ofile << "SCALING VELOCITY: " << vel_sa1 << endl;
            //jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1.z);

            //Adat[atom1].Ptp1 = Adat[atom1].Pt + vel_sa1;
            //Adat[atom2].Ptp1 = Adat[atom2].Pt - vel_sa1;
            vel_scales[atom1] += vel_sa1;
            //vel_scales[atom2] -= vel_sa1;

            //}
            //optfile->ofile << "Atom: " << i << " Ptm1: " << Adat[i].Ptm1 << " Pt: " << Adat[i].Pt << " Ptp1: " << Adat[i].Ptp1 << " VelCorr: " << vel_sa1 << " rN: " << rN << " Lambda: " << lambda <<"\n";
            //optfile->ofile << "Atom: " << atom2 << " Ptm1: " << Adat[atom2].Ptm1 << " Pt: " << Adat[atom2].Pt << " Ptp1: " << Adat[atom2].Ptp1 <<"\n\n";
        }
    }

    if (step % 500 == 0)
    {
        zeroLinearMomentum(data_mem,optfile,vel_scales);
    }

    for (int i = 0; i < N; ++i)
    {
        Adat[i].Ptp1 = Adat[i].Pt + vel_scales[i];
    }
};

/*____________________________________________________________________________
                   ----Linear Momentum Correction----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 1/21/2014  ----
                   ----Modified By:              ----
*/
void md_funcs::zeroLinearMomentum(MemHandler *data_mem,dataOutput* optfile,vector<jsm::vec3<double>> &vel)
{
    optfile->ofile << "\n|-------Scaling Bond Momentum-------|\n";
    if (!data_mem->coupled)
    {
        optfile->ofile << "Correcting Bond Momentums (Uncoupled)..." << "\n";
        for (int i = 0; i < K; ++i)
        {
            jsm::vec3<double> aP; //Calculate average linear momentum of the bond

            int atom1 = Bdat[i].atom1;
            int atom2 = Bdat[i].atom2;

            jsm::vec3<double> wv1 = vel[atom1];
            jsm::vec3<double> wv2 = vel[atom2];

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

            vel[atom1] = iVi1;
            vel[atom2] = iVi2;

            optfile->ofile << "\n|------------------------------|\n";
        }
    }
    else
    {
        optfile->ofile << "Correcting Molecular Momentum (Coupled)..." << "\n";
        jsm::vec3<double> aP; //Calculate average linear momentum of the molecule
        for (int i = 0; i < N; ++i)
        {
            aP.z += vel[i].z;
        }

        aP.z = aP.z/(double)N;//Average Momentum
        optfile->ofile << "Total Momentum Vector (Before Correction): " << aP << "\n";

        for (int i = 0; i < N; ++i)
        {
            vel[i].z = vel[i].z - aP.z;
        }
    }

    for (int i = 0; i < N; ++i)
    {
        double delta = 1.0E-10;
        if (abs(vel[i].x) < delta)
        {
            vel[i].x = 0.0f;
        }
        if (abs(vel[i].y) < delta)
        {
            vel[i].y = 0.0f;
        }
        if (abs(vel[i].z) < delta)
        {
            vel[i].z = 0.0f;
        }

    }
};

