/*
 * CSCues.h
 *
 *  Created on: Feb 2, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef CSCUES_H_
#define CSCUES_H_

#include "Cues.h"

class CSCues: public Cues {
public:
	CSCues(double div);
	virtual ~CSCues();
	void CSCuesImage(InputArray src, OutputArray dst, int featureType);
	void CSMapCalculation(InputArray src, OutputArray dst);
};

#endif /* CSCUES_H_ */
