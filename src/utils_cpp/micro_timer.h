#ifndef MICROTIMER_HPP
#define MICROTIMER_HPP

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                               Timer Class
//                  Holds timer variables and class functions
//      *************************************************************     //
class microTimer
{
    //--------------------------
    //Private Class Declarations
    //--------------------------

    double start_wall_time; //Holds start time
    long int wt_count; //Holds end time
    double accumtime;

    clock_t start_clock_time;
    long int ct_count;
    double accumclock;

    //------------------------------
    //Private Member Class Functions
    //------------------------------

    //Intakes a value of time in seconds and returns a string formmatted as:
    //Computation Wall Time: d days h hours m minutes s seconds
    string mk_time_string(double time_val);

public:
    microTimer()
    {
        Init();
    };

    //-----------------------------
    //Public Member Class Functions
    //-----------------------------
    // Initialie the class
    void Init(void);

    //Sets the start time for the timer
    void start_point(void);

    //Sets the end time for the timer
    void end_point(void);

    //Prints run_timer = end_time - start_time
    void print_wall_time(string message,ofstream &file1);

    //Prints the clock time
    void print_clock_time (string message,ofstream &file1);

    // Generic Time Print
    void print_generic(string message,ofstream &file1);

    //Resets the timer if needed
    void reset(void);
};

#endif

