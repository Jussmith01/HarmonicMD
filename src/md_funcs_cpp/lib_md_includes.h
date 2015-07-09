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
struct AtomicDataStore
{
    double AM; // Atomic mass
    double Kc; // Bond stiffness
    double r0; // Equilibrium bond distance

    std::vector<int> Nbonds; // Number of bonds on atom
    std::vector<int> Batom1; // Index of atomic bonds 1
    std::vector<int> Batom2; // Index of atomic bonds 2

    std::vector<double> Ptm1; // Position t-1
    std::vector<double> Pt; // Position t
    std::vector<double> Ptp1; // Position t+1

    std::vector<double> Vtm1; // Velocity t-1
    std::vector<double> Vt; // Velocity t
    std::vector<double> Vtp1; // Velocity t+1

    std::vector<double> Ftm1; // Force t-1
    std::vector<double> Ft; // Force t
    std::vector<double> Ftp1; // Force t+1
};

struct md_funcs
{
    //-------------------------
    //Struct Declarations
    //-------------------------
    //vector<mdAtomData> Adat;
    //vector<mdBondData> Bdat;
    std::vector<AtomicDataStore> adStore;

    //-------------------------
    //	  Data Storage
    //-------------------------

    //Initialized Variables
    float Ev_t;
    float Ev_tp1;

    float Tavg;

    //int *bond_vec; //bonding pattern MOVED TO mdBondData
    //int *bonds; //bonds per atom MOVED TO mdAtomData

    int N; // Number of atoms
    int K;
    float M;
    float dt;
    int steps;
    int step;

    float rmsdSum;
    int rmsdCount;
    vector<float> RMSD;

    //Constructor - Some basic initialization
    md_funcs()
    {
        Ev_t = 0;
        Tavg = 100;
        rmsdCount = 0;
        rmsdSum = 0;
    };

    //-------------------------
    //Struct Functions
    //-------------------------
    md_funcs(MemHandler *data_mem,dataOutput* optfile)
    {

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

