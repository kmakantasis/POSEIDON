/*
 * FrequencyFeatures.h
 *
 *  Created on: Feb 2, 2013
 *      Author: kostas makantasis
 */
#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef FREQUENCYFEATURES_H_
#define FREQUENCYFEATURES_H_

class FrequencyFeatures {
private:
	Mat imageGray;
public:
	FrequencyFeatures();
	virtual ~FrequencyFeatures();
	void FindFrequency(InputArray src, OutputArray dst);
};

#endif /* FREQUENCYFEATURES_H_ */
