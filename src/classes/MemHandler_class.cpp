#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../utils_cpp/lib_includes.h"
#include "class_header.h"

// ********************************************************************* //
// ****************cuMemHandler CLASS MEMBER FUNCTIONS****************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Allocates Necessary Memory ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 11/25/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the program run.
*/
void MemHandler::mem_alloc (dataOutput* optfile)
{
    N = ipt_parms.num_atoms;

    atomsMemSize = N * sizeof(atoms);
    memory_req = N * sizeof(int) + sizeof(double) + 11 * N * sizeof(double) + atomsMemSize;//

    //atom_data = new atoms [N];
    pos_vec = new double [N];
    vlc_vec = new double [N];
    atom_vec = new int [N];
}

/*____________________________________________________________________________
                   ----Calculate The Total Mass    ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 11/25/2014   ----
                   ----Modified By:                ----
*/
void MemHandler::Calc_Total_Mass()
{
    Mtot = 0;
    for (int i = 0; i < N; ++i)
    {
        Mtot += atom_data[i].AtomMass();
    }
    cout << "Mtot: " << Mtot << "\n";
}

/*____________________________________________________________________________
                   ----Input CM Shift		   ----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 11/25/2014   ----
                   ----Modified By:                ----
*/
void MemHandler::Input_CM_Shift (dataOutput* optfile)
{
    Calc_Total_Mass();//Set Mtot class variable

    double num,CM;
    num = 0.00;

    for (int i = 0; i < N; ++i)
    {
        //for (int j = 0; j < 3; ++j)
        //{
            num += atom_data[i].AtomMass() * atom_data[i].pos_xyz[2];
        //}
    }

    //for (int j = 0; j < 3; ++j)
    //{
        CM = num / (double) Mtot;
    //}

    //optfile->ofile << "N: " << N << "\n";
    //optfile->ofile << "K: " << K << "\n";

    for (int i = 0; i < N; ++i)
    {
        //for (int j = 0; j < 3; ++j)
        //{
            //optfile->ofile << "i: " << i << " j: " << j << "\n";
            atom_data[i].pos_xyz[2] = atom_data[i].pos_xyz[2] - CM;
        //}
    }
}

/*____________________________________________________________________________
                   ----Allocates Necessary Memory ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 11/25/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the program run.
*/
void MemHandler::mem_alloc_bonds (dataOutput* optfile)
{
    K = ipt_parms.num_bonds;

    bonddata = new int [2 * K];
    k = new double [K];
    r0 = new double [K];

    memory_req += 4 * K * sizeof(int) + 2 * K * sizeof(double);//

    //Set Unit Size for Memory Required output
    string unit = "B";
    double div = 1;
    if (memory_req >= 1024 && memory_req < 1024 * 1024)
    {
        div = 1024;
        unit = "KB";
    }
    else if (memory_req >= 1024 * 1024 && memory_req < 1024 * 1024 *1024)
    {
        div = 1024 * 1024;
        unit = "MB";
    }
    else if (memory_req >= 1024 * 1024 * 1024)
    {
        div = 1024 * 1024 * 1024;
        unit = "GB";
    }

    optfile->ofile << "Memory Required for Computation: " << memory_req / (double)div << unit <<"\n";
}

/*____________________________________________________________________________
                   ----Frees All Memory Allocations ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 11/25/2014     ----
                   ----Modified By:                 ----
Frees any memory allocated by this class.
*/
void MemHandler::mem_free ()
{
    delete [] bonddata;
    bonddata = NULL;
    delete [] atom_vec;
    atom_vec = NULL;
    delete [] pos_vec;
    pos_vec = NULL;
    delete [] vlc_vec;
    vlc_vec = NULL;
    delete [] k;
    k = NULL;
    delete [] r0;
    r0 = NULL;
    //delete [] atom_data; atom_data = NULL;
}

/*____________________________________________________________________________
                   ----Move Data to Working Memory----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 11/25/2014   ----
                   ----Modified By:               ----
*/
void MemHandler::move_wk_mem ()
{
    for (int i = 0; i < N; ++i)
    {
        atom_vec[i] = (int)atom_data[i].atomic_num;
        //for (int j = 0; j < 3; ++j)
        //{
            pos_vec[i] = atom_data[i].pos_xyz[2];
            vlc_vec[i] = atom_data[i].vlc_xyz[2];
        //}
    }
};

/*____________________________________________________________________________
                   ----DefineSystems		     ----
                   ----Author: Justin Smith          ----
                   ----Date Modified: 2/25/2015      ----
                   ----Modified By:                  ----
*/
void MemHandler::DefineSystems(dataOutput* optfile)
{
    int Nsys = systems.size();
    bool cpld = systems[0].coupled;
    coupled = cpld;

    int NumAtoms=0;
    int NumBonds=0;

    if(cpld)//if coupled
    {
        for (int i = 0; i < Nsys; ++i)
        {
            NumBonds += systems[i].numOscillators;
        }

        NumAtoms = NumBonds + 1;

        ipt_parms.num_atoms = NumAtoms;
        ipt_parms.num_bonds = NumBonds;

        mem_alloc(optfile);
        mem_alloc_bonds(optfile);

        DefineCoupled(optfile);

    }
    else if(!cpld)  //if not coupled
    {
        for (int i = 0; i < Nsys; ++i)
        {
            NumAtoms += 2 * systems[i].numOscillators;
            NumBonds += systems[i].numOscillators;
        }

        ipt_parms.num_atoms = NumAtoms;
        ipt_parms.num_bonds = NumBonds;

        mem_alloc(optfile);
        mem_alloc_bonds(optfile);

        DefineUnCoupled(optfile);
    }
};

/*____________________________________________________________________________
                   ----Define Coupled System         ----
                   ----Author: Justin Smith          ----
                   ----Date Modified: 2/25/2015      ----
                   ----Modified By:                  ----
*/
void MemHandler::DefineCoupled(dataOutput* optfile)
{
    int Nsys = systems.size();

    float pos = 0;
    int atomcnt = 0;

    //Define Atoms
    for (int i = 0; i < Nsys; ++i)
    {
        int NumAtoms = systems[i].numOscillators;

        if (i==0)
        {
            ++NumAtoms;
        };

        for (int j=0; j<NumAtoms; ++j)
        {
            atoms atom;

            atom.atomic_num = systems[i].AN;

            float r0 = systems[i].r0;

            atom.pos_xyz[0] = 0.0f;
            atom.pos_xyz[1] = 0.0f;
            atom.pos_xyz[2] = pos + r0;

            atom.vlc_xyz[0] = 0.0f;
            atom.vlc_xyz[1] = 0.0f;
            atom.vlc_xyz[2] = 0.0f;

            pos += r0;

            atom_data.push_back(atom);
            //cout << "SAVING ATOM: " << atomcnt << endl;
            ++atomcnt;
        }
    }

    Input_CM_Shift(optfile);
    move_wk_mem();

    //Define Bonds
    int bondcount = 0;
    for (int i = 0; i < Nsys; ++i)
    {
        int NumBonds = systems[i].numOscillators;

        for (int j=0; j<NumBonds; ++j)
        {
            //float r0 = systems[i].r0;
            bonddata[bondcount] = bondcount;
            bonddata[K + bondcount] = bondcount + 1;

            this->k[bondcount] = systems[i].k;
            this->r0[bondcount] = systems[i].r0;

            ++bondcount;
        }
    }
};

/*____________________________________________________________________________
                   ----Define Uncoupled System       ----
                   ----Author: Justin Smith          ----
                   ----Date Modified: 2/25/2015      ----
                   ----Modified By:                  ----
*/
void MemHandler::DefineUnCoupled(dataOutput* optfile)
{
    int Nsys = systems.size();

    //float pos = 0;
    int atomcnt = 0;

    //Define Atoms
    for (int i = 0; i < Nsys; ++i)
    {
        int poscnt = 0;
        int NumAtoms = 2 * systems[i].numOscillators;

        for (int j=0; j<NumAtoms; ++j)
        {
            atoms atom;

            atom.atomic_num = systems[i].AN;

            float r0 = systems[i].r0;

            atom.pos_xyz[0] = (float)i;
            atom.pos_xyz[1] = (float)poscnt;
            atom.pos_xyz[2] = (float)(j % 2) * r0;

            atom.vlc_xyz[0] = 0.0f;
            atom.vlc_xyz[1] = 0.0f;
            atom.vlc_xyz[2] = 0.0f;

            atom_data.push_back(atom);
            ++atomcnt;

            if (j % 2 == 1)
            {
                ++poscnt;
            }
        }
    }

    Input_CM_Shift(optfile);
    move_wk_mem();

    //Define Bonds
    int bondcount = 0;
    for (int i = 0; i < Nsys; ++i)
    {
        int NumBonds = systems[i].numOscillators;

        for (int j=0; j<NumBonds; ++j)
        {
            //float r0 = systems[i].r0;
            bonddata[bondcount] = 2 * bondcount;
            bonddata[K + bondcount] = 2 * bondcount + 1;

            this->k[bondcount] = (double)systems[i].k;
            this->r0[bondcount] = (double)systems[i].r0;

            ++bondcount;
        }
    }
};

