/*
 * ColorFeatures.h
 *
 *  Created on: Feb 4, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef COLORFEATURES_H_
#define COLORFEATURES_H_

class ColorFeatures {
private:
	Mat imageCIELAB;
public:
	ColorFeatures();
	virtual ~ColorFeatures();
	void FindColor(InputArray src, OutputArray dst, int featreType);
};

#endif /* COLORFEATURES_H_ */
