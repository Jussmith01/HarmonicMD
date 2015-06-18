#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>

//**Included custom headers**
#include "lib_includes.h"
#include "../classes/class_header.h"

using namespace std;

// ********************************************************************* //
// ****************DATA OUTPUT CLASS MEMBER FUNCTIONS******************* //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Create and Open Output File----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
This opens the output file by the name of 'filename'. If the output filename
was not included as the second argument upon command execution then the
output will be sent to a file named DEFAULT_OUTPUT.opt.
*/
void dataOutput::set_output (char *filename)
{
    //Define output name string
    string OutputFile;

    //Check if name is not included as the argument, then set output name.
    if (filename == NULL)
    {
        OutputFile = "DEFAULT_OUTPUT.opt";
    }
    else
    {
        OutputFile = filename;
    }

    //Open the output file by given filename.
    ofile.open(OutputFile.c_str());
}

/*____________________________________________________________________________
                   ----Close Output File          ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
This saves the data and closes the output file. Do this before ending the
program. When input int close_type = 0 closes normally, when 1, abnormal
termination.
*/
void dataOutput::close_output (int close_type)
{
    switch(close_type)
    {
    case 0:
    {
        ofile << "Normal Termination\n";
        break;
    }
    case 1:
    {
        ofile << "Abnormal Termination\n";
        break;
    }
    }

    ofile.close();
    CloseGraphOutputs();
}

/*____________________________________________________________________________
                   ----Prints Initial Data to Output----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 7/21/2014     ----
                   ----Modified By:                 ----
This prints the initial data to the output file.
*/
void dataOutput::initial_data_printer (MemHandler data_mem)
{
    ofile << "\n|------------------------------------------------------|\n";
    ofile << "|---------------------Start Atom Data------------------|\n";
    ofile << "|------------------------------------------------------|\n";

    ofile << "Molecular Data Loaded -- Number of Atoms = " << data_mem.ipt_parms.num_atoms << "\n\n";
    ofile << "Input Position and Velocity Data:\n";
    for (int i = 0; i < data_mem.ipt_parms.num_atoms; ++i)
    {
        ofile << "  Atom(" << i << ") ";
        ofile << " " << data_mem.atom_data[i].AtomLetter().c_str() << " ";
        if (data_mem.atom_data[i].AtomLetter().compare("error") == 0)
        {
            ostringstream error;
            error << "Error - Atomic number " << data_mem.atom_data[i].atomic_num << " is unsupported by this software.";
            Prog_err_quit(error.str(),data_mem);
        }
        ofile << " x(" << data_mem.atom_data[i].pos_xyz[0] << ") ";
        ofile << "y(" << data_mem.atom_data[i].pos_xyz[1] << ") ";
        ofile << "z(" << data_mem.atom_data[i].pos_xyz[2] << ") ";

        data_mem.atom_data[i].CartToSphere();

        //ofile << "Spherical Coords: R(" << data_mem.atom_data[i].pos_sph[0] << ") ";
        //ofile << "phi(" << data_mem.atom_data[i].pos_sph[1] << ") ";
        //ofile << "theta(" << data_mem.atom_data[i].pos_sph[2] << ") \n";

        ofile << "Vx(" << data_mem.atom_data[i].vlc_xyz[0] << ") ";
        ofile << "Vy(" << data_mem.atom_data[i].vlc_xyz[1] << ") ";
        ofile << "Vz(" << data_mem.atom_data[i].vlc_xyz[2] << ")\n";
    }

    ofile << "\nInput Bonding Data:\n";
    int L = data_mem.ipt_parms.num_bonds;
    for (int i = 0; i < L; ++i)
    {
        ofile << "  " << data_mem.bonddata[i] << "<->" << data_mem.bonddata[i + L] <<"\n";
    }
    ofile << "|------------------------------------------------------|\n";
    ofile << "|----------------------End Atom Data-------------------|\n";
    ofile << "|------------------------------------------------------|\n";
}

/*____________________________________________________________________________
                   ----Error Printer and Ender      ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 7/21/2014     ----
                   ----Modified By:                 ----
This prints a closing error message and ends the program.
*/
void dataOutput::Prog_err_quit(string Message,MemHandler data_mem)
{
    ofile << "\n" << Message << "\n";
    ofile << "Freeing Memory Allocations...\n";
    data_mem.mem_free();
    close_output(1);
    exit(1);
}

/*____________________________________________________________________________
                   ----Prepare Other Output Files ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 6/9/2015    ----
                   ----Modified By:               ----
This class prepares all other output files.
*/
void dataOutput::PrepareGraphOutputs(string dir)
{
    //ofstream graph[7];
    stringstream os[8];

    /*os[0] << dir << "Etot_graph.dat";
    os[1] << dir << "Vtot_graph.dat";
    os[2] << dir << "Ktot_graph.dat";
    os[3] << dir << "Ttot_graph.dat";
    os[4] << dir << "AvgRMSD.dat";
    os[5] << dir << "InstRMSD.dat";
    os[6] << dir << "mdout.xyz";*/

    os[0] << "Etot_graph.dat";
    os[1] << "Vtot_graph.dat";
    os[2] << "Ktot_graph.dat";
    os[3] << "Ttot_graph.dat";
    os[4] << "AvgRMSD.dat";
    os[5] << "InstRMSD.dat";
    os[6] << "mdout.xyz";
    os[7] << "Ztrans_graph.dat";

    for (int i = 0; i < 8; ++i)
    {
        graph[i].open(os[i].str().c_str());
    }
}

/*____________________________________________________________________________
                   ----Prepare Other Output Files ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 6/9/2015    ----
                   ----Modified By:               ----
This class prepares all other output files.
*/
void dataOutput::CloseGraphOutputs()
{
    for (int i = 0; i < 8; ++i)
    {
        graph[i].close();
    }
}
