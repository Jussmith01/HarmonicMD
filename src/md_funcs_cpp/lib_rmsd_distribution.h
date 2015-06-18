#ifndef dist_header
#define dist_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include "../math/jmath.h"

using namespace std;
// ********************************************************************* //
// **********************DEFINE MD WORKING CLASSES********************** //
// ********************************************************************* //
class rmsd_distribution
{
    string filename;
    vector<float> RMSD;

    double sig;//standard deviation
    double mu;//mean
    int N;

public:
    rmsd_distribution(vector<float> RMSD,string dir)
    {
        filename = "Prmsdvsrmsd.dat";
        this->RMSD = RMSD;
        N = RMSD.size();

        cout << "Calculating RMSD Distribution...\n";
        calculate_mu();
        calculate_sig();
        print_graph();
    };

    void calculate_mu()
    {
        double RMSDtot=0;
        for (int i=0; i<N; ++i)
        {
            RMSDtot += RMSD[i];
        }

        mu = RMSDtot / (double)N;
        cout << "Mean: " << mu << endl;
    };

    void calculate_sig()
    {
        double vTot=0;
        for (int i=0; i<N; ++i)
        {
            vTot += (RMSD[i] - mu) * (RMSD[i] - mu);
        }

        sig = sqrt(vTot  / (double)N);
        cout << "Std. Dev.: " << sig << endl;
    }

    void print_graph()
    {
        ofstream graph;

        int NUMPTS = 2000;
        float MIN = 0.0;
        float MAX = 2.5;

        double INC = (MAX - MIN) / (double)NUMPTS;

        double C = 1 / (sig * 2.50662827463);

        graph.open(filename.c_str());

        for (int i=0; i<N; ++i)
        {
            double x = i * INC;
            graph << x << "   " << C * exp(-((x-mu)*(x-mu))/(2.0f * sig * sig)) << endl;
        }

        graph.close();
    };
};

#endif

