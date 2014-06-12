/*
 * RightAngleFeatures.h
 *
 *  Created on: Feb 4, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>
#include "EdgeFeatures.h"

using namespace cv;

#ifndef RIGHTANGLEFEATURES_H_
#define RIGHTANGLEFEATURES_H_

class RightAngleFeatures {
private:
	Mat imageEdges;
public:
	RightAngleFeatures();
	virtual ~RightAngleFeatures();
	void FindRightAngles(InputArray src, OutputArray dst);
};

#endif /* RIGHTANGLEFEATURES_H_ */
