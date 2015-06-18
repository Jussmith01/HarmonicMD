#ifndef main_header
#define main_header

//**C/C++ Library Headers**
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <time.h>
#include <vector>

struct System
{
	bool coupled;
	int numOscillators;
	int AN;
	float k;//Spring Constant
	float r0;//Equilibrium 

	System(std::string coupled) 
	{
		if (coupled.compare("COUPLED")==0)
			{this->coupled=true;}
		else
			{this->coupled=false;}
	};

	System() {};
};

#endif

