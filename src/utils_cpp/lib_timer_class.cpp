#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <time.h>
#include "lib_includes.h"

using namespace std;

void timer::set_timer ()
{
    time(&start_time);
    t = clock();
}

void timer::end_timer ()
{
    time(&end_time);
    run_time = difftime(end_time,start_time);
    t = clock() - t;
    CLOCK_TIME = t / (double)CLOCKS_PER_SEC;
}

void timer::print_time (string message,ofstream &file1)
{
    file1 << message.c_str() << mk_time_string((double)run_time).c_str() << "\n";
}

void timer::print_clock_time (string message,ofstream &file1)
{
    file1 << message.c_str() << mk_time_string((double)CLOCK_TIME).c_str() << "\n";
}

void timer::print_time_s (string message,ofstream &file1)
{
    file1 << message.c_str() << run_time << "s\n";
}

void timer::reset ()
{
    start_time = 0;
    end_time = 0;
    run_time = 0;
    t = 0;
    CLOCK_TIME = 0;
}

string timer::mk_time_string(double time_val)
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
}

