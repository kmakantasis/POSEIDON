/*
 * CvANNMLPExt.h
 *
 *  Created on: Mar 10, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

#ifndef CVANNMLPEXT_H_
#define CVANNMLPEXT_H_

class CvANN_MLP_Ext: public CvANN_MLP {
public:
	CvANN_MLP_Ext();
	virtual ~CvANN_MLP_Ext();
};

#endif /* CVANNMLPEXT_H_ */
