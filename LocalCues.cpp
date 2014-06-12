/*
 * LocalCues.cpp
 *
 *  Created on: Jan 30, 2013
 *      Author: kostas makantasis
 */

#include "LocalCues.h"

LocalCues::LocalCues(double div) {this->blockDiv = div;}

LocalCues::~LocalCues() {}

void LocalCues::LocalCuesImage(InputArray src, OutputArray dst, int featureType){

	Mat srcMat = src.getMat();
	original = srcMat.clone();

	PyramidFeature(src, imagePyr_1, featureType);
	BlockDivision(imagePyr_1, imagePyr_1);

	CreatePyramid(src, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_2, featureType);
	BlockDivision(imagePyr_2, imagePyr_2);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_3, featureType);
	BlockDivision(imagePyr_3, imagePyr_3);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_4, featureType);
	BlockDivision(imagePyr_4, imagePyr_4);

	PyramidSummation(dst);
}
