#ifndef md_includes_header
#define md_includes_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include "../classes/class_header.h"
#include "../math/jmath.h"

using namespace std;
// ********************************************************************* //
// **********************DEFINE MD WORKING CLASSES********************** //
// ********************************************************************* //
struct mdAtomData
{
    double AM; // Atomic Mass

    jsm::vec3<double> Pi;

    jsm::vec3<double> Ptm1;
    jsm::vec3<double> Pt;
    jsm::vec3<double> Ptp1;

    jsm::vec3<double> Vtm1;
    jsm::vec3<double> Vt;
    jsm::vec3<double> Vtp1;

    int NBonds; // Number of atoms bonded
    vector<int> BA; // Bonded atoms
    vector<double> Kc; // Bonded atoms
    vector<double> r0; // Bonded atoms

    jsm::vec3<double> Ftm1;
    jsm::vec3<double> Ft;
    jsm::vec3<double> Ftp1;

    mdAtomData(float AM)
    {
        this->AM = (double)AM;
        NBonds = 0;
    }
};

struct mdBondData
{
    double Kc;
    double r0;

    int atom1; //Defines the atoms involved in the bond
    int atom2;

    mdBondData(float Kc)
    {
        this->Kc = (double)Kc;
    }
};

class md_funcs
{
    public:
    //-------------------------
    //Struct Declarations
    //-------------------------
    vector<mdAtomData> Adat;
    vector<mdBondData> Bdat;

    //Initialized Variables
    float Ev_t;
    float Ev_tp1;

    float Tavg;

    //int *bond_vec; //bonding pattern MOVED TO mdBondData
    //int *bonds; //bonds per atom MOVED TO mdAtomData

    int N;
    int K;
    float M;
    float dt;
    int steps;
    int step;

    bool heated;
    bool pisaved;

    float rmsdSum;
    int rmsdCount;
    vector<float> RMSD;

    //Constructor - Some basic initialization
    md_funcs()
    {
        std::cout << "CONSTRUCTOR!\n";
        Ev_t = 0;
        Tavg = 100;
        rmsdCount = 0;
        rmsdSum = 0;
        heated=false;
        pisaved=false;
    };

    //-------------------------
    //Struct Functions
    //-------------------------
    md_funcs(MemHandler *data_mem,dataOutput* optfile)
        :md_funcs()
    {
        //ofstream graph[7];
        //stringstream os[7];

        //os[0] << data_mem->ipt_parms.data_dir << "Etot_graph.dat";
        //os[1] << data_mem->ipt_parms.data_dir << "Vtot_graph.dat";
        //os[2] << data_mem->ipt_parms.data_dir << "Ktot_graph.dat";
        //os[3] << data_mem->ipt_parms.data_dir << "Ttot_graph.dat";
        //os[4] << data_mem->ipt_parms.data_dir << "AvgRMSD.dat";
        //os[5] << data_mem->ipt_parms.data_dir << "InstRMSD.dat";
        //os[6] << data_mem->ipt_parms.data_dir << "mdout.xyz";

        //for (int i = 0; i < 7; ++i)
        //{
            //graph[i].open(os[i].str().c_str());
            //graph[i].close();
        //}
    };

    ~md_funcs() {};

    void md_mem_alloc(MemHandler *data_mem,dataOutput* optfile);
    //void Calculate_Initial_Contraction (dataOutput* optfile);
    //void md_mem_free (dataOutput* optfile);
    void mv_starting_vectors(MemHandler *data_mem,dataOutput* optfile);

    void CM_Calc(MemHandler *data_mem,dataOutput* optfile);
    //float rad_r3v(float *pos_vec,int atom1, int atom2,int sqroot);
    //double travel_dist(int atom);
    void calc_bonds_per_atom(dataOutput* optfile);
    void calc_forces(MemHandler *data_mem,dataOutput* optfile);
    void verlet_integration(MemHandler *data_mem,dataOutput* optfile);
    void shift_pos_vec(dataOutput* optfile);
    //void print_bond_dists(dataOutput* optfile);
    double calc_E_total(MemHandler *data_mem,dataOutput* optfile);
    double calc_potential(int bond,dataOutput* optfile);
    double calc_kinetic(int atom,dataOutput* optfile);
    void produce_md_out(MemHandler *data_mem,dataOutput* optfile);
    void scale_velocities(MemHandler *data_mem,dataOutput* optfile,double T);
    void zeroLinearMomentum(MemHandler *data_mem,dataOutput* optfile,vector< jsm::vec3<double> > &vel);
    void Velocity_Initialization(MemHandler *data_mem,dataOutput* optfile);

    //RMSD
    void calculate_rmsd(MemHandler *data_mem,dataOutput* optfile);
    void ShiftCM(vector< jsm::vec3<double> > &position,dataOutput* optfile);
};

// ********************************************************************* //
// **********************DEFINE EXTERNAL FUNCTIONS********************** //
// ********************************************************************* //

//
string Harmonic_MD_main (MemHandler *data_mem,dataOutput *optfile);

#endif

