#ifndef mem_includes
#define mem_includes
//C/C++ includes
#include "../mainHeader.h"

//Custom includes
#include "../utils_cpp/lib_includes.h"

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //

//________________________________________________________________________//
//      *************************************************************     //
//                      Program Memory Handler Class
//      This class hold all needed memory handling declarations and
//	actions, such as memory allocations, transfers, to and from 
//	device and host, in the form of class member functions.
//      *************************************************************     //
class MemHandler
{
	public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
		//Parameters handler
		parameters ipt_parms;
	
		//Molecular data handler
		size_t atomsMemSize;
		vector<atoms> atom_data;

		//Working Memory Storage
		long int memory_req;
		double *pos_vec; //Position Vector -- Linearly stored as x1,y1,z1,x2,y2,z2,...,xn,yn,zn
		double *vlc_vec; //Velocity Vector -- Linearly stored as x1,y1,z1,x2,y2,z2,...,xn,yn,zn
		double *k;//Spring Constants
		double *r0;//Equilibrium bond distances
		int *bonddata;//Holds the bonding data
		int *atom_vec;//vector holding atomic numbers
		float Mtot;
		int N,K;//N = Number of atoms, K = Number of bonds, Mtot = Mass total

		bool coupled;
		vector<System> systems;
	
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
		//Allocated Memory
		void mem_alloc(dataOutput* optfile);
		void mem_alloc_bonds(dataOutput* optfile);

		//Free All Memory
		void mem_free(void);

		//Move working memory
		void move_wk_mem(void);

		//Build Systems
		void DefineSystems(dataOutput* optfile);
		void DefineCoupled(dataOutput* optfile);
		void DefineUnCoupled(dataOutput* optfile);

		//Shifts the input Coordinates
		void Input_CM_Shift(dataOutput* optfile);
		void Calc_Total_Mass(void);
};

// ********************************************************************* //
// ***************************DEFINE FUNCTIONS**************************** //
// ********************************************************************* //

#endif
