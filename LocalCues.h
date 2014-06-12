/*
 * LocalCues.h
 *
 *  Created on: Jan 30, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef LOCALCUES_H_
#define LOCALCUES_H_

#include "Cues.h"

class LocalCues: public Cues {
public:
	LocalCues(double div);
	virtual ~LocalCues();
	void LocalCuesImage(InputArray src, OutputArray dst, int featureType);
};

#endif /* LOCALCUES_H_ */
