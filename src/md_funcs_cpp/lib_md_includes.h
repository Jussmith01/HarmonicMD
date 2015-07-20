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
// ******************************************************************** //
struct mdData
{
    // Atomic Data
    std::vector< double > AM;

    std::vector< jsm::vec3<double> > Pi;

    std::vector< jsm::vec3<double> > Ptm1;
    std::vector< jsm::vec3<double> > Pt;
    std::vector< jsm::vec3<double> > Ptp1;

    std::vector< jsm::vec3<double> > Vtm1;
    std::vector< jsm::vec3<double> > Vt;
    std::vector< jsm::vec3<double> > Vtp1;

    vector<int> Nbonds;
    std::vector< std::vector<int> > BA;
    std::vector< std::vector<double> > aKc;
    std::vector< std::vector<double> > ar0;

    std::vector< jsm::vec3<double> > Ftm1;
    std::vector< jsm::vec3<double> > Ft;
    std::vector< jsm::vec3<double> > Ftp1;

    // Bond Data
    std::vector< double > bKc;
    std::vector< double > br0;

    std::vector< unsigned int > atom1;
    std::vector< unsigned int > atom2;

    // Clear all allocated memory
    void clear()
    {
        if (!AM.empty()) {AM.clear();}

        if (!Pi.empty()) {Pi.clear();}
        if (!Ptm1.empty()) {Ptm1.clear();}
        if (!Pt.empty()) {Pt.clear();}
        if (!Ptp1.empty()) {Ptp1.clear();}

        if (!Vtm1.empty()) {Vtm1.clear();}
        if (!Vt.empty()) {Vt.clear();}
        if (!Vtp1.empty()) {Vtp1.clear();}

        if (!Nbonds.empty()) {Nbonds.clear();}

        if (!BA.empty()) {BA.clear();}
        if (!aKc.empty()) {aKc.clear();}
        if (!ar0.empty()) {ar0.clear();}

        if (!Ftm1.empty()) {Ftm1.clear();}
        if (!Ft.empty()) {Ft.clear();}
        if (!Ftp1.empty()) {Ftp1.clear();}

        if (!bKc.empty()) {bKc.clear();}
        if (!br0.empty()) {br0.clear();}

        if (!atom1.empty()) {atom1.clear();}
        if (!atom2.empty()) {atom2.clear();}
    };
};

class md_funcs
{
    //-------------------------
    //Struct Declarations
    //-------------------------
    mdData dataStore;

public:

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

    int SaveSteps;
    int mdSaveSteps;

    bool heated;
    bool pisaved;

    float rmsdSum;
    int rmsdCount;
    vector<float> RMSD;

    //Constructor - Some basic initialization
    md_funcs()
    {
        Ev_t = 0;
        Tavg = 100;
        SaveSteps=5;
        mdSaveSteps=50;
        rmsdCount = 0;
        rmsdSum = 0;
        heated=false;
        pisaved=false;
    };

    //-------------------------
    //Struct Functions
    //-------------------------
    md_funcs(MemHandler *data_mem,dataOutput* optfile) : md_funcs()
    {

    };

    ~md_funcs() {dataStore.clear();};

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
    double calc_potential(int bond);
    double calc_kinetic(int atom,double ifact);
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

