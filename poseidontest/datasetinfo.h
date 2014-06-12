// Class to hold information about the datasets in the test-set folder.
// Paris Kaimakis 29 Jan 2013.


#ifndef _DATASETINFO_H_
#define _DATASETINFO_H_

#include <iostream>


#define NUM_DATASETS	10


#define F_QUANTIFY




struct s_dataset
{
	std::string szFilename;
	std::string szDescription;
};
typedef struct s_dataset t_dataset;



struct s_info
{
	t_dataset aDataset[NUM_DATASETS];
};
typedef struct s_info t_info;






#endif // _DATASETINFO_H_
