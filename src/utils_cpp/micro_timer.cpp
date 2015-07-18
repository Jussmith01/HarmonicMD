#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string.h>
#include <string>
#include <time.h>
#include <omp.h>
#include "micro_timer.h"

using namespace std;

void microTimer::Init ()
{
    wt_count=0;
    ct_count=0;
    accumtime=0;
    accumclock=0;
};

void microTimer::start_point ()
{
    start_wall_time = omp_get_wtime();
    start_clock_time = clock();
};

void microTimer::end_point ()
{
    accumtime+=omp_get_wtime()-start_wall_time;
    accumclock+=clock()-start_clock_time / (double)CLOCKS_PER_SEC;

    ++wt_count;
    ++ct_count;

    start_clock_time=0;
    start_wall_time=0;
};

void microTimer::print_wall_time (string message,ofstream &file1)
{
    //file1 << message.c_str() << mk_time_string((double)accumtime/(double)wt_count).c_str() << "\n";
    file1 << message.c_str() << mk_time_string((double)accumtime).c_str() << "\n";
    reset();
};

void microTimer::print_clock_time (string message,ofstream &file1)
{
    file1 << message.c_str() << mk_time_string((double)accumclock/(double)ct_count).c_str() << "\n";
    reset();
};

void microTimer::print_generic(string message,ofstream &file1)
{
    stringstream ss1;
    ss1 << message << " Wall Time: ";

    //stringstream ss2;
    //ss2 << message << " Clock Time: ";

    print_wall_time(ss1.str(),file1);
    //print_clock_time(ss2.str(),file1);
};

void microTimer::reset ()
{
    wt_count=0;
    ct_count=0;

    accumtime=0;
    accumclock=0;
};

string microTimer::mk_time_string(double time_val)
{
    int days, hours, minutes;
    double seconds;

    days = floor(time_val/86400);
    hours = floor((time_val-(days*86400))/3600);
    minutes = floor((time_val - (days * 86400)-(hours * 3600))/60);
    seconds = (double)(time_val - (days * 86400) - (hours * 3600) - (minutes * 60));

    ostringstream time_str;
    time_str << days << " days " << hours  << " hours " << minutes << " minutes " << seconds  << " seconds";
    return time_str.str();
};

