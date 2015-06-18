//**C/C++ included libraries**
#include "mainHeader.h"

//**Included CUDA libraries**

//**Included custom headers**
#include "utils_cpp/lib_includes.h"
#include "classes/class_header.h"
#include "md_funcs_cpp/lib_md_includes.h"

using namespace std;

int main (int argc, char *argv[])
{
    //********************************************//
    //               Declarations                 //
    //********************************************//
    char *inputfile = argv[1];
    char *opfilename = argv[2];
    MemHandler data_mem;//Define Memory Handler
    dataOutput optfile(opfilename);//Define Output Handler
    timer wac_timer;//Define Timer Handler

    //********************************************//
    //              Preamble Programs             //
    //********************************************//
    wac_timer.set_timer();//Set the program wall timer

    //********************************************//
    //              Read Inputs                   //
    //********************************************//
    read_input(inputfile,&data_mem,&optfile);//Read the input file

    //********************************************//
    //            Execution Premble               //
    //********************************************//
    data_mem.ipt_parms.define_parms(&optfile);
    optfile.PrepareGraphOutputs(data_mem.ipt_parms.data_dir);
    data_mem.DefineSystems(&optfile);
    optfile.initial_data_printer(data_mem);//Print input data to output file

    //********************************************//
    //          Main Program Execution            //
    //********************************************//
    string err_str = Harmonic_MD_main(&data_mem,&optfile);

    //********************************************//
    //              Post Run Programs             //
    //********************************************//
    wac_timer.end_timer();//End the program wall timer
    wac_timer.print_clock_time("Clock Time: ",optfile.ofile);//Print the wall time to output
    wac_timer.print_time("Wall Time: ",optfile.ofile);//Print the wall time to output
    optfile.close_output(0);//Close the output file
    data_mem.mem_free();//Free all memory

    return 0;
}
