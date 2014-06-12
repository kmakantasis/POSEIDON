/*
 * GlobalCues.h
 *
 *  Created on: Jan 31, 2013
 *      Author: kostas
 */

#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef GLOBALCUES_H_
#define GLOBALCUES_H_

#include "Cues.h"

class GlobalCues: public Cues {
public:
	GlobalCues(double div);
	virtual ~GlobalCues();
	void GlobalCuesImage(InputArray src, OutputArray dst, int featureType);
	void GlobalMapCalculation(InputArray src, OutputArray dst);
};

#endif /* GLOBALCUES_H_ */
