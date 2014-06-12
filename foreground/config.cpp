// Config reader/handler implementation
// Paris Kaimakis 29 Jan 2013


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "config.h"

using std::cout;
using std::endl;



double string_to_double( const std::string& s )
 {
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
 } 







void readConfigFile(const std::string szPath, t_settings* pSets)
{
	cout << "\tReading config file... " ;
	cout.flush();


	std::string szFilename = szPath + "config.txt";
	std::ifstream fin(szFilename.c_str(), std::ios::in);

	//cout << szFilename << endl;
	if(!fin.good())
		cout << "Unable to open config file." << endl;
	


	struct s_cfgVar
	{
		std::string szDescription;
		int iType;
		double var;
	};
	typedef struct s_cfgVar t_cfgVar;

	t_cfgVar aVar[NUM_CONFIG_VARS];
	aVar[0].szDescription = "framewidth";
	aVar[0].iType = E_INT;
	aVar[1].szDescription = "frameheight";
	aVar[1].iType = E_INT;
	aVar[2].szDescription = "minvariance";
	aVar[2].iType = E_DOUBLE;
	aVar[3].szDescription = "minclustersize";
	aVar[3].iType = E_DOUBLE;
	aVar[4].szDescription = "DeltaStDev";
	aVar[4].iType = E_DOUBLE;
	aVar[5].szDescription = "thSigma";
	aVar[5].iType = E_DOUBLE;
	aVar[6].szDescription = "sigmaMin";
	aVar[6].iType = E_DOUBLE;


	while(!fin.eof())
	{
		char czLine[200];
		
		fin.getline(czLine, 100);
		std::string szLine = czLine;
		//cout << szLine<< endl;


		
		int pos = szLine.find('#');
		if(pos!=std::string::npos)
			szLine.erase(pos, szLine.length());
		//cout << szLine << endl;
		


		for(int iVar=0; iVar<NUM_CONFIG_VARS; ++iVar)
		{
			std::string szDescription = aVar[iVar].szDescription;
			double& var = aVar[iVar].var;

			int pos = szLine.find(szDescription );
			if(pos!=std::string::npos)
			{
				// Get the next word, convert it into number.
				szLine.erase(0,szDescription.length());
				
				
				switch(aVar[iVar].iType)
				{
				case E_INT:
					var = atoi(szLine.c_str());
					break;
				case E_DOUBLE:
					var = string_to_double(szLine);
					break;
				}


				//cout << "iV=" << iVar << " var=" << var << endl;

			}
		}
		//*/

	}

	fin.close();



	pSets->U				= (int) aVar[0].var;
	pSets->V				= (int) aVar[1].var;
	pSets->minVariance		= aVar[2].var;
	pSets->minClusterSize	= aVar[3].var;
	pSets->DeltaStDev		= aVar[4].var;
	pSets->thSigma			= aVar[5].var;
	pSets->sigmaMin			= aVar[6].var;

	cout << "Done!" << endl;
}
