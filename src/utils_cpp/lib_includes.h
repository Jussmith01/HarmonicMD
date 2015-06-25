#ifndef class_header
#define class_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //
//Predeclarations to prevent co-dependency compilation errors
class atoms;
class parameters;
class output_file_set;
class MemHandler;

//________________________________________________________________________//
//      *************************************************************     //
//				Data Output Class
//      Holds everything needed to output the data too the output files.
//      *************************************************************     //
class dataOutput
{
    //-------------------------
    //Public Class Declarations
    //-------------------------
public:
    dataOutput (char *filename)
    {
        set_output(filename);
    };

    // Declare the output stream:
    ofstream ofile;
    int verbose;

    ofstream graph[8];
    //-----------------------------
    //Public Member Class Functions
    //-----------------------------

    //Creates/Opens the output file:
    void set_output(char *filename);

    //Closes the output file:
    void close_output(int close_type);

    //Prints the initial input data to the output:
    void initial_data_printer(MemHandler data_mem); //Prints the initial input data from the input

    //Program to help with error handling printing. Will cut down on lines of code later.
    void Prog_err_quit(string Message,MemHandler data_mem);

    // Setup other outputs
    void PrepareGraphOutputs(string dir);

    // Setup other outputs
    void CloseGraphOutputs();
};

//________________________________________________________________________//
//      *************************************************************     //
//                              Atomic Data Class
//      Holds the data for atoms along with conversion and transfor-
//	mation functions. An array of these give molecular data.
//      *************************************************************     //
class atoms
{

    //-------------------------
    //Public Class Declarations
    //-------------------------
public:
    atoms (const atoms &instance)
    {
        this->atomic_num = instance.atomic_num;
        this->pos_xyz[0] = instance.pos_xyz[0];
        this->pos_xyz[1] = instance.pos_xyz[1];
        this->pos_xyz[2] = instance.pos_xyz[2];
        this->vlc_xyz[0] = instance.vlc_xyz[0];
        this->vlc_xyz[1] = instance.vlc_xyz[1];
        this->vlc_xyz[2] = instance.vlc_xyz[2];
    };

    atoms () {};

    //Holds the atomic number:
    int atomic_num;

    //Holds the cartesian position vector of the atom:
    float pos_xyz[3]; // 0 = x; 1 = y; 2 = z; 3 = k;

    //Holds the spherical position vector of the atom:
    float pos_sph[3]; //0 = r; 1 = theta; 2 = phi

    //Holds the velocity vector of the atom:
    float vlc_xyz[3]; // 0 = Vx; 1 = Vy; 2 = Vz;

    //-----------------------------
    //Public Member Class Functions
    //-----------------------------

    //Converts the cartesian vector in pos_xyz to a spherical vector in pos_sph:
    void CartToSphere(void);

    //Converts the sperical vector in pos_sph to a cartesian vector in pos_xyz:
    void SphereToCart(void);

    //Returns the periodic element letter code for the given value in atomic_num:
    string AtomLetter(void);

    //Returns the mass of the atom:
    double AtomMass(void);
};

//________________________________________________________________________//
//      *************************************************************     //
//                            Parameter Data Class
//      	    Holds any parameters from the input file.
//      *************************************************************     //
class parameters
{
    //-------------------------
    //Public Class Declarations
    //-------------------------
public:
    //This array holds the input parameters, specificed below:
    /*
      ---------------------------------------
      |Array Index   |      Description	    |
      ---------------------------------------
      | parmchars[0] |MD Temperature	    |
      |     [1]      |Verbosity		        |
      |     [2]      |Time Step		        |
      |	    [3]      |Steps			        |
      |	    [4]      |Data Directory	    |
      ---------------------------------------
    */
    string parmchars[7];

    //Holds the total number of atoms in the molecule:
    int num_atoms;
    int num_bonds;
    int steps;
    int heatSteps; //Steps for heating the system
    double temp;
    double dt;
    double tau; //Brendensen Thermostat tau value

    string data_dir;

    //Output verbosity
    int opt_verb;

    //
    void define_parms(dataOutput *optfile);
    //void produce_out_dir(dataOutput *optfile);
};

//________________________________________________________________________//
//      *************************************************************     //
//                               Timer Class
//                  Holds timer variables and class functions
//      *************************************************************     //
class timer
{
    //--------------------------
    //Private Class Declarations
    //--------------------------

    time_t start_time; //Holds start time
    time_t end_time; //Holds end time
    double run_time; //Holds run time = end_time - start_time
    clock_t t;
    double CLOCK_TIME;

    //------------------------------
    //Private Member Class Functions
    //------------------------------

    //Intakes a value of time in seconds and returns a string formmatted as:
    //Computation Wall Time: d days h hours m minutes s seconds
    string mk_time_string(double time_val);

public:
    //-----------------------------
    //Public Member Class Functions
    //-----------------------------

    //Sets the start time for the timer
    void set_timer(void);

    //Sets the end time for the timer
    void end_timer(void);

    //Prints run_timer = end_time - start_time
    void print_time(string message,ofstream &file1);

    //Prints the clock time
    void print_clock_time (string message,ofstream &file1);

    //Prints run_timer = end_time - start_time
    void print_time_s(string message,ofstream &file1);

    //Resets the timer if needed
    void reset(void);
};
// ********************************************************************* //
// **********************DEFINE EXTERNAL FUNCTIONS********************** //
// ********************************************************************* //

//This function reads the input file. Call with:
// read_input([Name of input file], &[parameters class declaration], &[a NULL atoms class declaration])
void read_input (char *inputfile,MemHandler *data_mem,dataOutput *optfile);

#endif

