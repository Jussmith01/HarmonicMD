//Included Libraries
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <cstring>
#include <cmath>

//Included Headers
#include "lib_includes.h"

using namespace std;

// ********************************************************************* //
// ****************ATOMS CLASS MEMBER FUNCTIONS************************* //
// ********************************************************************* //

/*____________________________________________________________________________
		----Cartesian to Spherical Converter----
		----Author: Justin Smith	    ----
                ----Date Modified: 7/21/2014        ----
                ----Modified By:                    ----
This atoms class member function changes the current cartesian coords loaded
into pos_xyz[i] array into spherical coords (degrees) in the pos_sph[i] array. 
*/
void atoms::CartToSphere ()
{

	float x,y,z,x2,y2,z2,tmp1,tmp2,tmp3;
	float degrees = 57.2957795;// <------- MAKE THIS A GLOBALLY DEFINED CONSTANT

	x = pos_xyz[0];	
	y = pos_xyz[1];
	z = pos_xyz[2];

	x2 = pow(x,2);
	y2 = pow(y,2);
	z2 = pow(z,2);
	
	tmp1 = sqrt(x2 + y2 + z2);// Calculate distance
	tmp2 = atan2(sqrt(x2 + y2),z); // Calculate theta
	tmp3 = atan2(y,x); // Calculate phi

        pos_sph[0] = tmp1; //Save calculated r from cartesian
        pos_sph[1] = tmp2 * degrees; //Save calculated theta from cartesian
        pos_sph[2] = tmp3 * degrees; //Save calculated phi from cartesian
}

/*____________________________________________________________________________
                ----Spherical To Cartesian Converter----
                ----Author: Justin Smith            ----
                ----Date Modified: 7/21/2014        ----
                ----Modified By:                    ----
This atoms class member function changes the current spherical coords (degrees) loaded
into pos_sph[i] array into cartesian coords in the pos_xyz[i] array. 
*/

void atoms::SphereToCart ()
{

        float r,theta,phi,tmp1,tmp2,tmp3;
        float degrees = 57.2957795;//<-------- MAKE THIS A GLOBALLY DEFINED CONSTANT

        r = pos_sph[0]; 
        theta = pos_sph[1] / degrees;
        phi = pos_sph[2] / degrees;

	tmp1 = r * sin(theta) * cos(phi); //Calculate x
	tmp2 = r * sin(theta) * sin(phi); //Calculate y
	tmp3 = r * cos(theta);	//Calculate z
		
        pos_xyz[0] = tmp1; //Save calculated x from spherical
       	pos_xyz[1] = tmp2; //Save calculated y from spherical
        pos_xyz[2] = tmp3; //Save calculated z from spherical
}

/*____________________________________________________________________________
                ----Atomic Letter Code Selector Program----
                ----Author: Justin Smith               ----
		----Date Modified: 7/21/2014	       ----
		----Modified By:		       ----
This atoms class member function selects the proper letter code for the given 
atomic number defined. The main purpose of this is simply for output, but 
also, any unsupported atoms should not be included here. This is so that the 
unsupported atom error is thrown and the program is stopped at this point with 
the correct error message.
*/	
string atoms::AtomLetter ()
{
	string Letter;
	switch (atomic_num)
	{
		case 1: {Letter = "H";break;}
		case 5: {Letter = "B";break;}
		case 6: {Letter = "C";break;}
		case 7: {Letter = "N";break;}
		case 8: {Letter = "O";break;}
		case 9: {Letter = "F";break;}
		case 10: {Letter = "Ne";break;}
		default: {Letter = "error";break;}
	}
	return Letter;
}

/*____________________________________________________________________________
                ----Atomic Letter Code Selector Program----
                ----Author: Justin Smith               ----
                ----Date Modified: 7/21/2014           ----
                ----Modified By:                       ----
This atoms class member function selects the proper letter code for the given 
atomic number defined. The main purpose of this is simply for output, but 
also, any unsupported atoms should not be included here. This is so that the 
unsupported atom error is thrown and the program is stopped at this point with 
the correct error message.
*/
double atoms::AtomMass()
{
       double mass;
        switch (atomic_num)
        {
                case 1: {mass = 1.6737236E-27;break;}
                case 5: {mass = 1.79520863E-26;break;}
                case 6: {mass = 1.99442348E-26;break;}
                case 7: {mass = 1.99442348E-26;break;}
                //case 7: {mass = 2.32587369E-26;break;}
                case 8: {mass = 2.65676264E-26;break;}
                case 9: {mass = 3.15475880E-26;break;}
                case 10: {mass = 3.35091773E-26;break;}
                default: {exit(0);break;}
        }
        return mass;
}

