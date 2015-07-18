#ifndef vTools_class_header
#define vTools_class_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <random>
#include <time.h>
#include <omp.h>
#include "../math/jmath.h"

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

namespace tools
{
template<typename T>
inline void precisionsetd(jsm::vec3<T>& vec,T delta)
{
    if (abs(vec[0]) < delta)
        vec[0] = 0.0f;

    if (abs(vec[1]) < delta)
        vec[1] = 0.0f;

    if (abs(vec[2]) < delta)
        vec[2] = 0.0f;
};

class RandomRealVal
{
    std::default_random_engine generator;
    std::vector<int> array;
    int index;

public:
    RandomRealVal() {};

    RandomRealVal(int w,int i) {Setup(w,i);};

    void Setup(int w,int i)
    {
        time_t Time;
        time(&Time);
        int seedOffset=(int)Time;

        array.resize(w);

        int t = (int)omp_get_wtime()+i;
        std::seed_seq seed = {seedOffset,t,i+100};
        seed.generate(array.begin(),array.end());//Seed the generator
        index = 0;
    };

    double GenRandReal(double high,double low)
    {
        generator.seed(array[index]);//Seed the generator
        std::uniform_real_distribution<double> distribution(low,high);//Setup the distribution
        double RN = (double)distribution(generator);//Denerate the random number
        ++index;//Increase seed offset
        return RN;
    };
};
}
// ********************************************************************* //
// **********************DEFINE EXTERNAL FUNCTIONS********************** //
// ********************************************************************* //

#endif

