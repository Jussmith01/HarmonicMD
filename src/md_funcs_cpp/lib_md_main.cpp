#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iostream>
#include <iomanip>

#include "../utils_cpp/lib_includes.h"
#include "../classes/class_header.h"
#include "lib_md_includes.h"
#include "lib_rmsd_distribution.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                              Main MD Function
//      Carrys out the MD steps.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 12/16/2014  ----
                   ----Modified By:               ----
*/

extern string Harmonic_MD_main(MemHandler *data_mem,dataOutput *optfile)
{
    optfile->ofile << "\n|------------------------------------------------------|\n";
    optfile->ofile << "|----------------Begin MD Calculations-----------------|\n";
    optfile->ofile << "|------------------------------------------------------|\n";

    //********************************************//
    //               Declarations                 //
    //********************************************//
    md_funcs md_tools(data_mem,optfile);

    //********************************************//
    //              Preamble Programs             //
    //********************************************//
    md_tools.md_mem_alloc(data_mem,optfile); //Allocates required memory for md cycles
    md_tools.mv_starting_vectors(data_mem,optfile); // populates pos_t and bond_vec
    //Need program to calculate an initial molecule contraction. This sets pos_vec(t-1)
    //to begin the MD simulation integration.

    //********************************************//
    //             Main MD Execution              //
    //********************************************//
    optfile->ofile << "\nBeginning MD Steps...\n";
    for (int i = 0; i < md_tools.steps; ++i)//i is step number
    {
        //cout << "STEP: " << i << "\n";
        md_tools.step = i;
        optfile->ofile << "|-----------------Step (" << i << ")-----------------|\n";
        //Next, need to calculate forces. Can begin with initial forces = 0 since at time
        //step t = t0 r = req therefore F = 0.

        //Scale the velocities to match temperature
        //md_tools.scale_velocities(data_mem,optfile);

        md_tools.calc_forces(data_mem,optfile); //Populates F_t

        //Next, need to integrate via Verlet.
        md_tools.verlet_integration(data_mem,optfile); //determines pos_tp1 via the Verlet algorithm.
        //Next, need a thermostat to maintain an Etot based on temp T.

        md_tools.CM_Calc(data_mem,optfile);

        if (i != 0)
        {
            double Etot = md_tools.calc_E_total(data_mem,optfile);

            if (i % 15 == 0)
            {
                md_tools.produce_md_out(data_mem,optfile);
            }
            //md_tools.produce_md_out(data_mem,optfile);
            md_tools.calculate_rmsd(data_mem,optfile);
            //md_tools.print_bond_dists(optfile);
        }

        md_tools.shift_pos_vec(optfile);
        optfile->ofile << "|--------------------------------------------|\n\n";
        //`cout << "\n";
    }

    rmsd_distribution calcDist(md_tools.RMSD,data_mem->ipt_parms.data_dir);

    string errstr = "NO";

    //********************************************//
    //              Post Run Programs             //
    //********************************************//

    optfile->ofile << "\n|------------------------------------------------------|\n";
    optfile->ofile << "|------------------End MD Calculations-----------------|\n";
    optfile->ofile << "|------------------------------------------------------|\n\n";

    return errstr;
}
