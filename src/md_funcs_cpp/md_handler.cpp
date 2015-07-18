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
#include "../utils_cpp/lib_various_tools.h"


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

    dataStore.Pi.resize(N);

    dataStore.Ftm1.resize(N);
    dataStore.Ft.resize(N);
    dataStore.Ftp1.resize(N);

    //Copy stating atom data
    dataStore.AM.resize(N);
    for(int i=0;i<N;++i)
    {
        dataStore.AM[i]=data_mem->atom_data[i].AtomMass();
    }

    //Produce bonding data stuff
    dataStore.bKc.resize(K);
    //memcpy(&dataStore.bKc[0],&data_mem->k[0],K*sizeof(double));

    tools::RandomRealVal RN(K,clock());

    std::ofstream rng;
    rng.open("rnGenPlot.dat");

    for (int i = 0 ; i < K; ++i)
    {
        dataStore.bKc[i]=RN.GenRandReal(100.0,10.0);
        rng << i << "   " << dataStore.bKc[i] << std::endl;
    }
    rng.close();

    dataStore.br0.resize(K);
    memcpy(&dataStore.br0[0],&data_mem->r0[0],K*sizeof(double));

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
    //Copy initial position and velocity vector from Memhandler
    dataStore.Pt.resize(N);
    memcpy(&dataStore.Pt[0],&data_mem->pos_vec[0],N*sizeof(jsm::vec3<double>));

    dataStore.Vt.resize(N);
    memcpy(&dataStore.Vt[0],&data_mem->vlc_vec[0],N*sizeof(jsm::vec3<double>));

    //Copy bonding data from Memhandler
    dataStore.atom1.resize(K);
    memcpy(&dataStore.atom1[0],&data_mem->bonddata[0],K*sizeof(int));

    dataStore.atom2.resize(K);
    memcpy(&dataStore.atom2[0],&data_mem->bonddata[K],K*sizeof(int));

    // Setup atomic bonding data
    calc_bonds_per_atom(optfile);

    // Check if velocities need to be initialized
    double tMag = 0.0;
    for (int i = 0; i < N; ++i)
        tMag += jsm::magnitude(dataStore.Vt[i]);

    // Initialize velocities if needed
    if (tMag < 1.0E-6)
    {
        optfile->ofile << "No starting velocities detected, running velocity initialization.\n";
        Velocity_Initialization(data_mem,optfile);
    }
    else
    {
        optfile->ofile << "Starting velocities detected, no need for velocity initialization.\n";
    }

    //Copy initial position and velocity vector from Memhandler
    dataStore.Ptm1.resize(N);
    dataStore.Ptp1.resize(N);
    for (int i = 0 ; i < N; ++i)
    {
        dataStore.Ptm1[i] = dataStore.Pt[i] - dataStore.Vt[i];
        dataStore.Ptp1[i] = dataStore.Pt[i] + dataStore.Vt[i];
    }

    //Print the position vectors
    optfile->ofile << "\nPositions: \n";
    for (int i=0;i<N;++i)
    {
        optfile->ofile << " tm1: " << i << dataStore.Ptm1[i] << "\n";
        optfile->ofile << " t: " << i << dataStore.Pt[i] << "\n";
        optfile->ofile << " tp1: " << i << dataStore.Ptp1[i] << "\n";
    }

    //Print the velocity vectors
    optfile->ofile << "\nInitial Velocities: \n";
    for (int i = 0; i < N; ++i)
    {
        optfile->ofile << " " << i << ") AM: " << dataStore.AM[i] << dataStore.Vt[i] << "\n";
    }

    optfile->ofile << "\nBonding: \n";
    for (int i = 0; i < K; ++i)
    {
        optfile->ofile << " " << i << ": ";
        optfile->ofile << dataStore.atom1[i] << " <-> ";
        optfile->ofile << dataStore.atom2[i] << "\n";
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

    //tools::RandomRealVal RN(K,clock());

    for (int i = 0; i < N; ++i)
    {
        optfile->ofile << "\n|----Calc velocities for atom: " << i << "----|\n";

        //jsm::vec3<double> vtemp(0.0,0.0,pow(-1,i) * Vscale);
        dataStore.Vt[i] = jsm::vec3<double>(0.0,0.0,pow(-1,i) * Vscale);

        optfile->ofile << "vT(" << i << "): " << dataStore.Vt[i] << endl;
    }

    optfile->ofile << "\n|-------Scaling Bond Momentum-------|\n";
    for (int i = 0; i < K; ++i)
    {
        jsm::vec3<double> aP; //Calculate average linear momentum of the bond

        int atom1 = dataStore.atom1[i];
        int atom2 = dataStore.atom2[i];

        jsm::vec3<double> wv1 = dataStore.Vt[atom1];
        jsm::vec3<double> wv2 = dataStore.Vt[atom2];

        wv1 = jsm::UniformScalar(wv1,dataStore.AM[atom1]);
        wv2 = jsm::UniformScalar(wv2,dataStore.AM[atom2]);

        aP = jsm::UniformScalar(wv1 + wv2,0.5);
        optfile->ofile << "Bond(" << i << "): Total Bond Momentum Vector (Before Correction): " << aP << "\n";

        jsm::vec3<double> iP1 = jsm::UniformScalar(dataStore.Vt[atom1],dataStore.AM[atom1]) - aP;
        jsm::vec3<double> iP2 = jsm::UniformScalar(dataStore.Vt[atom2],dataStore.AM[atom2]) - aP;

        dataStore.Vt[atom1] = UniformScalar(iP1, 1/(double)dataStore.AM[atom1]);
        dataStore.Vt[atom2] = UniformScalar(iP2, 1/(double)dataStore.AM[atom2]);
    }

    for (int i = 0; i < N; ++i)
        tools::precisionsetd(dataStore.Vt[i],1.0E-14);
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
            num[j] += dataStore.AM[i] * dataStore.Pt[i][j];
        }
    }

    double invM = 1.0/(double)M;
    for (int j = 0; j < 3; ++j)
    {
        CM[j] = num[j] * invM;
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
    dataStore.BA.resize(N);
    dataStore.aKc.resize(N);
    dataStore.ar0.resize(N);

    dataStore.Nbonds.resize(N);
    memset(&dataStore.Nbonds[0],0,N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < K; ++j)
        {
            if (i == (int)dataStore.atom1[j])
            {
                dataStore.BA[i].push_back( dataStore.atom2[j] );
                dataStore.aKc[i].push_back( dataStore.bKc[j] );
                dataStore.ar0[i].push_back( dataStore.br0[j] );
                ++dataStore.Nbonds[i];
            }

            if (i == (int)dataStore.atom2[j])
            {
                dataStore.BA[i].push_back( dataStore.atom1[j] );
                dataStore.aKc[i].push_back( dataStore.bKc[j] );
                dataStore.ar0[i].push_back( dataStore.br0[j] );
                ++dataStore.Nbonds[i];
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
    //#pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        //Determine Total Force Components (for each F[x] = F[0],F[y] = F[1],F[z] = F[2])
        jsm::vec3<double> Ftot;

        int atomA = i; // Atom A
        jsm::vec3<double> pA = dataStore.Pt[atomA];

        //Calculate forces for each atom
        for (int l = 0; l < dataStore.Nbonds[i]; ++l)
        {
            int atomB = dataStore.BA[i][l];
            double Kc = dataStore.aKc[atomA][l];
            double r0 = dataStore.ar0[atomA][l];
            jsm::vec3<double> pB = dataStore.Pt[atomB];

            jsm::vec3<double> x = pA - pB;
            jsm::vec3<double> x0 = jsm::normalize(x) * r0;

            jsm::vec3<double> rv = x - x0;

            //cout << "TEST\n";

            //Zero values below precision
            tools::precisionsetd(rv,1.0E-14);

            //optfile->ofile << "Kc: " << Kc << " x: " << x << " x0: " << x0 << " x-x0: " << rv << "\n";

            jsm::vec3<double> Fvec = (rv) * (-Kc);
            Ftot += Fvec;
        }

        dataStore.Ft[i] = Ftot;
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
        double AM = (double)dataStore.AM[i];

        jsm::vec3<double> av = dataStore.Ft[i] * (dt2/(double)AM);
        jsm::vec3<double> t1 = dataStore.Pt[i] * 2.0;

        tools::precisionsetd(av,1.0E-14);

        double x = dataStore.Pt[i][0];
        double y = dataStore.Pt[i][1];
        double z = t1[2] - dataStore.Ptm1[i][2] + av[2];
        dataStore.Ptp1[i] = jsm::vec3< double >(x,y,z);
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
        dataStore.Ptm1[i] = dataStore.Pt[i];
        dataStore.Pt[i] = dataStore.Ptp1[i];
    }
}

/*____________________________________________________________________________
                   ----Calculate Potential       ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
double md_funcs::calc_potential(int bond)
{
    int atom1 = dataStore.atom1[bond];
    int atom2 = dataStore.atom2[bond];

    double r_T = jsm::magnitude(dataStore.Pt[atom1] - dataStore.Pt[atom2]);

    float Kc = dataStore.bKc[bond];
    float r0 = dataStore.br0[bond];

    return Kc * (r_T - r0) * (r_T - r0);
}

/*____________________________________________________________________________
                   ----Calculate Kinetic         ----
                   ----Author: Justin Smith      ----
                   ----Date Modified: 12/17/2014 ----
                   ----Modified By:              ----
*/
double md_funcs::calc_kinetic(int atom,double ifact)
{
    double Tdist = jsm::magnitude(dataStore.Ptm1[atom] - dataStore.Ptp1[atom]);
    double AM = (double)dataStore.AM[atom];

    double v = Tdist * ifact;

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
    double Vtot = 0.0;

    //#pragma omp parallel for default(shared) reduction(+:Vtot)
    for (int i = 0; i < K; ++i)
    {
        Vtot = calc_potential(i);
    }
    Vtot = 0.5 * Vtot;

    //Calculate the total kinetic energy of the atomic motions
    double Ktot = 0.0;
    double ifact = 1.0 / (double)(2.0 * dt);
    //#pragma omp parallel for default(shared) firstprivate(ifact) reduction(+:Ktot)
    for (int i = 0; i < N; ++i)
    {
        Ktot += calc_kinetic(i,ifact);
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
    if (heated)
    {
        if (step % 50 == 0)
        {
            int atom1 =

            dataStore.atom1[0];
            int atom2 = dataStore.atom2[0];

            float radius = jsm::magnitude(dataStore.Pt[atom1] - dataStore.Pt[atom2]);

            //cout << "SHOULD BE SAVING DATA!\n";

            optfile->graph[0] << radius << "  " << Etot << "\n";
            optfile->graph[1] << radius << "  " << Vtot << "\n";
            optfile->graph[2] << radius << "  " << Ktot << "\n";
            optfile->graph[3] << step * dt << "  " << T << "\n";
        }
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
    if (heated)
    {
        optfile->graph[6] << N <<"\n\n";

        for (int i = 0; i < N; ++i)
        {
            optfile->graph[6] << data_mem->atom_data[i].AtomLetter() << "    " << dataStore.Pt[i][0] << "      " << dataStore.Pt[i][1] << "      " << dataStore.Pt[i][2] << "\n";
        }
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
    if (heated)
    {
        vector<jsm::vec3<double>> position;
        position.resize(N);
        ShiftCM(position,optfile);

        ++rmsdCount;
        float rmsd = 0;
        for (int i = 0; i < N; ++i)
        {
            float d = position[i][2] - dataStore.Pi[i][2];
            rmsd += d * d;
        }

        rmsd = sqrt(rmsd / (float)N);
        rmsdSum += rmsd;

        //Save data after heating for statistics
        RMSD.push_back(rmsd);

        //Save RSMD Graphs
        if(step % 50 == 0)
        {
            //std::cout << "Saving RMSD Data...\n";
            optfile->graph[4] << step * dt * 1.0E12 << "  " << rmsdSum / (float)rmsdCount << "\n";
            optfile->graph[5] << step * dt * 1.0E12 << "  " << rmsd << "\n";
        }
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
        Ztotal += dataStore.Pt[i][2] * dataStore.AM[i];
    }
    double Zcm = Ztotal / (double)M;

    optfile->ofile << "Zcm: " << Zcm << endl;

    for (int i = 0; i < N; ++i)
    {
        position[i][2] = dataStore.Pt[i][2] - Zcm;
        //optfile->ofile << "Initial Pos: " << Adat[i].Pt << " Shifted Pos: " << position[i] << endl;

        if (step % 500 == 0)
        {
            dataStore.Ptm1[i][2] = dataStore.Ptm1[i][2] - Zcm;
            dataStore.Pt[i][2] = dataStore.Pt[i][2] - Zcm;
            dataStore.Ptp1[i][2] = dataStore.Ptp1[i][2] - Zcm;
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

    //std::default_random_engine generator;
    //std::uniform_real_distribution<double> distribution(-0.00001,0.00001);
    tools::RandomRealVal RN1(N,clock());

    //std::default_random_engine generator2;
    //std::uniform_real_distribution<double> distribution2(1.0,0.0);
    //tools::RandomRealVal RN2(K,clock());


    if (step < HS)
    {
        Tr = ceil(Tr * (step / (double)HS));

    }
    else
    {
        //std::cout << "HEATED!!\n";
        heated=true;
        if(!pisaved)
        {
            //Copy initial position  and velocity vector from Memhandler
            for (int i = 0 ; i < N; ++i)
            {
                    dataStore.Pi[i] = dataStore.Pt[i];
            }
            pisaved=true;
        }
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
            //if (distribution(generator2) < 0.5)
            //if (RN2.GenRandReal(1.0,0.0) < 0.5)
            //{
                //double rN = 0;
                //double rN = RN1.GenRandReal(0.00001,-0.00001);
                //int atom1 = i;
                int atom1 = dataStore.atom1[i];
                int atom2 = dataStore.atom2[i];

                jsm::vec3<double> vel_s1 = jsm::zScalar(dataStore.Ptp1[atom1] - dataStore.Ptm1[atom1],lambda/2.0f);

                jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1[2]);
                //jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1[2]+rN);
                //optfile->ofile << "SCALING VELOCITY: " << vel_sa1 << endl;
                //jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1.z);

                //Adat[atom1].Ptp1 = Adat[atom1].Pt + vel_sa1;
                //Adat[atom2].Ptp1 = Adat[atom2].Pt - vel_sa1;
                vel_scales[atom1] += vel_sa1;
                vel_scales[atom2] -= vel_sa1;

            //}
        }
    }
    else
    {
        double sumRandV=0;
        vector<double> rN;

        for (int i=0; i<N; ++i)
        {
            rN.push_back(RN1.GenRandReal(0.00001,-0.00001));
            sumRandV += rN[i];
        }

        sumRandV /= (double)N;
        for (int i=0; i<N; ++i)
        {
            rN[i] = rN[i] - sumRandV;
        }

        for (int i = 0; i < N; ++i)
        {
            int atom1 = i;
            jsm::vec3<double> vel_s1 = jsm::zScalar(dataStore.Ptp1[atom1] - dataStore.Ptm1[atom1],lambda/2.0f);

            jsm::vec3<double> vel_sa1(0.0,0.0,vel_s1[2]+rN[i]);
            vel_scales[atom1] += vel_sa1;
        }
    }

    if (step % 500 == 0)
    {
        zeroLinearMomentum(data_mem,optfile,vel_scales);
    }

    for (int i = 0; i < N; ++i)
    {
        dataStore.Ptp1[i] = dataStore.Pt[i] + vel_scales[i];
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

            int atom1 = dataStore.atom1[i];
            int atom2 = dataStore.atom2[i];

            jsm::vec3<double> swv1 = jsm::UniformScalar(vel[atom1],dataStore.AM[atom1]);
            jsm::vec3<double> swv2 = jsm::UniformScalar(vel[atom2],dataStore.AM[atom1]);

            aP += swv1 + swv2;
            //aP += swv2;

            jsm::vec3<double> aPa = jsm::UniformScalar(aP,0.5);
            optfile->ofile << "Bond(" << i << "): Total Bond Momentum Vector (Before Correction): " << aPa << "\n";

            jsm::vec3<double> iP1 = jsm::UniformScalar(dataStore.Vt[atom1],dataStore.AM[atom1]);
            jsm::vec3<double> iP2 = jsm::UniformScalar(dataStore.Vt[atom2],dataStore.AM[atom2]);

            jsm::vec3<double> iPi1 = iP1 - aPa;
            jsm::vec3<double> iPi2 = iP2 - aPa;

            jsm::vec3<double> iVi1 = UniformScalar(iPi1, 1/(double)dataStore.AM[atom1]);
            jsm::vec3<double> iVi2 = UniformScalar(iPi2, 1/(double)dataStore.AM[atom2]);

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
            aP[2] += vel[i][2];

        aP[2] = aP[2]/(double)N;//Average Momentum
        optfile->ofile << "Total Momentum Vector (Before Correction): " << aP << "\n";

        for (int i = 0; i < N; ++i)
            vel[i][2] = vel[i][2] - aP[2];
    }

    for (int i = 0; i < N; ++i)
        tools::precisionsetd(vel[i],1.0E-14);
};

