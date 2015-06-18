#ifndef vTools_class_header
#define vTools_class_header

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
class v_tools
{
        //--------------------------
        //Private Class Declarations
        //--------------------------

        //------------------------------
        //Private Member Class Functions
        //------------------------------

        public:
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
	string rm_whitespace(string line)
        {
		string s = line;
                size_t spot = s.find_first_of("\t\n");
                while (spot!=string::npos)
                {
                        //cout << spot << endl;
                        s.erase(spot);
                        spot = s.find_first_of("\t\n");
                }
                //cout << "string: " << s << endl;
                return s;
        };


};
// ********************************************************************* //
// **********************DEFINE EXTERNAL FUNCTIONS********************** //
// ********************************************************************* //

#endif

