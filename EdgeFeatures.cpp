/*
 * EdgeFeatures.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: kostas makantasis
 */

#include "EdgeFeatures.h"

EdgeFeatures::EdgeFeatures() {}

EdgeFeatures::~EdgeFeatures() {}

void EdgeFeatures::FindEdges(InputArray src, OutputArray dst, int minThres, int ratio, int window) {
	int maxThres = ratio*minThres;

	cvtColor(src, imageGray, CV_RGB2GRAY);
	Canny(imageGray, edgeCanny, minThres, maxThres, window);
	GaussianBlur(imageGray, imageGray, Size(3,3), 0, 0, BORDER_DEFAULT);
	Sobel(imageGray, edgeGx, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	convertScaleAbs(edgeGx, edgeGx);
	Sobel(imageGray, edgeGy, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
	convertScaleAbs(edgeGy, edgeGy);
	convertScaleAbs(edgeCanny, edgeCanny);
	addWeighted(edgeGx, 0.5, edgeGy, 0.5, 0, edgeGxGy);
	threshold(edgeCanny, edgeCanny, 128, 1, THRESH_BINARY);
	multiply(edgeGxGy, edgeCanny, dst);
}
