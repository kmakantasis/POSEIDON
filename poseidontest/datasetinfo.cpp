// Class to hold information about the datasets in the test-set folder.
// Paris Kaimakis 29 Jan 2013.



#include <iostream>
#include <sstream>
#include <fstream>
#include "datasetinfo.h"

using std::cout;
using std::endl;




void createDataSetInfo(t_info* pInfo)
{
	// Create, store info.
	
	for(int i=0; i<NUM_DATASETS; ++i)
	{
		t_dataset* pData = &pInfo->aDataset[i];
		
		std::stringstream ssNum;
		ssNum << i;
		pData->szFilename = ssNum.str() + ".avi";

	}



	t_dataset* aData = pInfo->aDataset;
	aData[0].szDescription = "Easy: Sunny morning, sailboats";
	aData[1].szDescription = "Easy: Sunny morning, sailboats, short, only 1 target appears";
	aData[2].szDescription = "Medium: Sunny morning, sailboats, very small targets also appear in background";
	aData[3].szDescription = "Easy: Sunny morning, big ship moving on the horizon";
	aData[4].szDescription = "Medium: Sunny morning, medium/small boat moving left-to-right, big ship moving right-to-left";
	aData[5].szDescription = "Medium/Hard: Sunny morning, sailboats, very small targets";
	aData[6].szDescription = "Hard: Average afternoon, boat";
	aData[7].szDescription = "Very Hard: Average afternoon, small boat, sometimes barely visible to a human eye";
	aData[8].szDescription = "Hard: Average afternoon, boat";
	aData[9].szDescription = "Hard: Changing lighting conditions, vandal";


}




void displayInfo(t_info* pInfo)
{
	for(int i=0; i<NUM_DATASETS; ++i)
	{
		std::stringstream ssNum ;
		ssNum << i;
		
		cout << ssNum.str() << ") " ;
		cout << pInfo->aDataset[i].szFilename << " - ";
		cout << pInfo->aDataset[i].szDescription << endl;
	}
	cout << endl;

}




void getPath(const std::string szPathWithExe, std::string& szPath)
{
	szPath = "";

	int pos = szPathWithExe.rfind('\\');
	for(int i=0; i<=pos; ++i)
	{
		szPath += szPathWithExe[i];
	}

	/*
#ifdef F_QUANTIFY
	//szPath = "C:\\Users\\Paris\\Documents\\Dev\\Poseidon\\poseidontest\\Release\\";
	szPath = "C:\\Users\\paris.kaimakis\\Documents\\Poseidon\\Code\\poseidontest\\Release\\";
#endif
	*/

	return;
}



