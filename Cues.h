/*
 * Cues.h
 *
 *  Created on: Jan 30, 2013
 *      Author: kostas makantasis
 */

#ifndef CUES_H_
#define CUES_H_

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;

class Cues {
protected:
	Mat original;
	Mat imagePyr_1;
	Mat imagePyr_2;
	Mat imagePyr_3;
	Mat imagePyr_4;
	Mat pyrTemp;
	Mat sumTemp, sumTempPyr12, sumTempPyr23, sumTempPyr34;
	Mat scaleTemp;
	double blockDiv;
public:
	Cues();
	virtual ~Cues();
	void BlockDivision(InputArray src, OutputArray dst);
	void CreatePyramid(InputArray src, OutputArray dst);
	int PyramidFeature(InputArray src, OutputArray dst, int featureType);
	void ColorMap(InputArray src, OutputArray dst);
	void PyramidSummation(OutputArray dst);
	//void CuesImage(InputArray src, int featureType);
};

#endif /* CUES_H_ */
