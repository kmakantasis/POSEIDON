/*
 * EdgeFeatures.h
 *
 *  Created on: Jan 29, 2013
 *      Author: kostas makantasis
 */
#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef EDGEFEATURES_H_
#define EDGEFEATURES_H_

class EdgeFeatures {
private:
	Mat imageGray;
	Mat edgeCanny;
	Mat edgeGx;
	Mat edgeGy;
	Mat edgeGxGy;
public:
	EdgeFeatures();
	virtual ~EdgeFeatures();
	void FindEdges(InputArray src, OutputArray Dst, int minThres, int ratio, int window);
};

#endif /* EDGEFEATURES_H_ */
