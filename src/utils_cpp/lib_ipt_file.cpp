/*******************************************************************
************File Designed to Read in XYZ Trajectory File************
********************************************************************
DESIGNED BY: JUSTIN S. SMITH
CONTACT E-MAIL: JSMITH48@CHEM.UFL.EDU

SAMPLE INPUT FILE:

-----------------STARTS BELOW THIS LINE---------------------
$PARAM
12              !par test int 1
4.1485          !par test flt 4
str_test        !par test string
$ENDPARAM
$COORD
        1       -1.12798434        1.545547           2.2132132
        6       -0.1224696362       -3.1174832504      -1.5015549412
        8        2.1012819183       4.0397895568      -1.1329875333
        6        3.2785162937       3.3434825848      -0.7159896664
        6       14.1137391171      14.5455261160      -3.6515096380
$ENDCOORD
$VELOC
        -1.3996393401   -8.0351111152   -4.4886276994
        7.5952796391    -2.4957179104     3.6110497172
        8.3281570155    -1.5098881501    -0.0410243452
        6.3406004079     3.2470414282    -7.2878335148
        4.4153733160     3.2574127691    -3.4879277197
$ENDVELOC
$END OF FILE$
-------------------ENDS ABOVE THIS LINE-----------------------
*********************************************************************/

//**C/C++ included libraries**
#include "../mainHeader.h"

//**Included custom headers**
#include "lib_includes.h"
#include "../classes/class_header.h"

using namespace std;

//*****************************************************************
//***********Obtains line number of file line delimiters***********
//*****************************************************************
unsigned int get_line_num (const char *search_string, string *INPUT_DATA, int NUM_LINES)
{
	unsigned int LINE_NUM=0,i=0;
	while (i < (unsigned int)NUM_LINES)
	{
		if (strcmp(search_string,INPUT_DATA[i].c_str()) == 0)
		{LINE_NUM = i;break;}
		i++;
	}
	return LINE_NUM;
}

//------------------------------------------
//      Parse A String by an : sign
//------------------------------------------
bool parseLine(string command,string &func,string &param)
{
        bool errchk = true;
        size_t pos = command.find(":");
        //cout << "POS: " << pos << endl;
        if (pos != string::npos)
        {
                func = command.substr(0,pos);
                param = command.substr(pos+1);
                errchk = false;
        } else {
                func = command;
                param = "NONE";
                errchk = false;
        }
        return errchk;
};

//*********************************************************************
//**********Save velocities of velocity vector to atoms class**********
//*********************************************************************
void save_system_setup (string *INPUT_DATA,MemHandler *data_mem,int LINE_BEG,int LINE_END)
{
        //int beg_sstr,end_sstr,first,last,sz_substr,i=0,error_chk;
        int i=0;
        string comp_str = "";
        for (int j = LINE_BEG + 1; j < LINE_END; ++j)
        {
		System system;
		string coupled;
		string wkcommand;

		//Set coupled
		if(parseLine(INPUT_DATA[j],coupled,wkcommand))
			{cout << "!!ERROR!!\n";};

		if(coupled.compare("UNCOUPLED") == 0)
			{system.coupled = false;}
		else if(coupled.compare("COUPLED") == 0)
			{system.coupled = true;}
		else
			{cout << "!!ERROR!! Unrecognised input: " << coupled << endl;}


		//Set number of oscillators
		string NumStr;
		if(parseLine(wkcommand,NumStr,wkcommand))
                        {cout << "!!ERROR!!\n";};

		system.numOscillators = atoi(NumStr.c_str());

                //Set number of Atomic Number
                if(parseLine(wkcommand,NumStr,wkcommand))
                        {cout << "!!ERROR!!\n";};

                system.AN = atoi(NumStr.c_str());

		//Set Spring const and equilibrium dist
                if(parseLine(wkcommand,NumStr,wkcommand))
                        {cout << "!!ERROR!!\n";};

                system.k = atof(NumStr.c_str());
                system.r0 = atof(wkcommand.c_str());

		cout << "SYSTEM DATA (" << i << "): " << coupled << " NO: " << system.numOscillators << " AN: " << system.AN << " k: " << system.k << " r0: " << system.r0 << endl;

		data_mem->systems.push_back(system);
		++i;
        }
}

//*******************************************************************
//**********************Read the input file**************************
//*******************************************************************
extern void read_input (char *inputfile,MemHandler *data_mem,dataOutput *optfile)
{

	//ipt_data *DATA;
	string line;
	ifstream iptfile (inputfile);
	int i=0,NUM_LINES=0;

	//COUNTS NUMBER OF LINES IN THE INPUT FILE FOR DYNAMIC MEMORY ALLOCATION
	if (iptfile.is_open())
        	{
                	while (!iptfile.eof())
                	{
				getline(iptfile,line);
				if (line == "$END OF FILE$") {break;}
				NUM_LINES++;
	                }
        iptfile.close();

        } else {optfile->ofile << "***error -- unable to open file***\n";exit(1);}

	//ALLOCATION OF TEMP DATA STRING ARRAY
	string *DATA = NULL;
	DATA = new string [NUM_LINES];

	//READ INPUT LINES INTO STRING ARRAY: DATA
	ifstream ipt2file (inputfile);
	if (ipt2file.is_open())
	{
		while (!ipt2file.eof())
		{
			getline(ipt2file,line);
			if (line != "$END OF FILE$")
			{
				DATA[i] = line;
			} else {break;}
			i++;
		}
	ipt2file.close();

	} else {optfile->ofile << "***error -- unable to open file***\n";exit(1);}
	//PARSE DATA FUNCTIONS FOLLOW

	unsigned int LINE_BEG, LINE_END,char_search=0;
	string test_string_beg,test_string_end;
	i = 0;
	while (i <= 3)
	{
		//cout << "TESTING ("<< i << ")!\n";
		switch (i)
		{
			case 0: {test_string_beg = "$PARAM";test_string_end = "$ENDPARAM";break;}
			case 1: {test_string_beg = "$SYSTEM";test_string_end = "$ENDSYSTEM";break;}
		}

                LINE_BEG = get_line_num(test_string_beg.c_str(),DATA,NUM_LINES);
                LINE_END = get_line_num(test_string_end.c_str(),DATA,NUM_LINES);

		switch (i)
		{
			case 0:
			{
				//cout << "PARMS\n";
				for (int j = (int)LINE_BEG + 1; j < (int)LINE_END; ++j)
				{
					//cout << "LINE(" << j << "): " << DATA[j] << "\n";
					char_search = DATA[j].find("!");
					data_mem->ipt_parms.parmchars[j-LINE_BEG-1] = DATA[j].substr(0,char_search);
					//cout << data_mem->ipt_parms.parmchars[j-LINE_BEG-1] << "\n";
				}
				break;
			}
                        case 1:
                        {
				//cout << "COORD\n";
                                //data_mem->ipt_parms.num_atoms = (LINE_END - LINE_BEG) - 1;
				//data_mem->mem_alloc(optfile);
				save_system_setup(DATA,data_mem,LINE_BEG,LINE_END);
				break;
                        }
		}
		i++;
	}

	//***********CLEAR FUNCTION MEMORY ALLOCATIONS***********
	delete [] DATA;
	DATA = NULL;
}
