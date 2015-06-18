#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <iostream>
#include "lib_includes.h"
#include "lib_various_tools.h"

using namespace std;

void parameters::define_parms(dataOutput *optfile)
{
	v_tools tools;
	optfile->ofile << "\n|------------------------------------------------------|\n";
	optfile->ofile << "|--------------------Input Parameters------------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n";
	istringstream(parmchars[0]) >> temp;
	optfile->ofile << "Temperature: " << temp << "\n";
	istringstream(parmchars[1]) >> opt_verb;
	optfile->ofile << "Output Verbosity: " << opt_verb << "\n";
	optfile->verbose = opt_verb;
        istringstream(parmchars[2]) >> dt;
        optfile->ofile << "Time Step: " << dt << "\n";
        istringstream(parmchars[3]) >> steps;
        optfile->ofile << "Total Steps: " << steps << "\n";
        data_dir = tools.rm_whitespace(parmchars[4]);
        optfile->ofile << "Data Directory: " << data_dir << "\n";
        istringstream(parmchars[5]) >> heatSteps;
        optfile->ofile << "Heating Steps: " << heatSteps << "\n";
        istringstream(parmchars[6]) >> tau;
        optfile->ofile << ": " << tau << "\n";	
	optfile->ofile << "|------------------------------------------------------|\n";
	optfile->ofile << "|-----------------End Input Parameters-----------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n\n";
}

