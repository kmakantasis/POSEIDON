/*
 * ColorFeatures.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: kostas makantasis
 */

#include "ColorFeatures.h"

ColorFeatures::ColorFeatures() {}

ColorFeatures::~ColorFeatures() {}

void ColorFeatures::FindColor(InputArray src, OutputArray dst, int featureType) {

	cvtColor(src, dst, CV_BGR2Lab);

	if(featureType == 41) {
		Scalar meanValuePyr_1 = mean(dst);
		absdiff(dst, meanValuePyr_1, dst);
	}

}

