/*
 * FrequencyFeatures.cpp
 *
 *  Created on: Feb 2, 2013
 *      Author: kostas makantasis
 */

#include "FrequencyFeatures.h"

FrequencyFeatures::FrequencyFeatures() {}

FrequencyFeatures::~FrequencyFeatures() {}

void FrequencyFeatures::FindFrequency(InputArray src, OutputArray dst) {

	Mat srcMat = src.getMat();
	GaussianBlur(srcMat, srcMat, Size(3,3), 0, 0, BORDER_DEFAULT );
	cvtColor(srcMat, imageGray, CV_RGB2GRAY);

	Laplacian(imageGray, dst, CV_16S, 3, 1, 0, BORDER_DEFAULT);
	convertScaleAbs(dst, dst);
}

