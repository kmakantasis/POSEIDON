/*
 * GlobalCues.cpp
 *
 *  Created on: Jan 31, 2013
 *      Author: kostas
 */

#include "GlobalCues.h"

GlobalCues::GlobalCues(double div) { this->blockDiv = div;}

GlobalCues::~GlobalCues() {}

void GlobalCues::GlobalCuesImage(InputArray src, OutputArray dst, int featureType) {

	Mat srcMat = src.getMat();
	original = srcMat.clone();

	PyramidFeature(src, imagePyr_1, featureType);
	BlockDivision(imagePyr_1, imagePyr_1);
	GlobalMapCalculation(imagePyr_1, imagePyr_1);

	CreatePyramid(src, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_2, featureType);
	BlockDivision(imagePyr_2, imagePyr_2);
	GlobalMapCalculation(imagePyr_2, imagePyr_2);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_3, featureType);
	BlockDivision(imagePyr_3, imagePyr_3);
	GlobalMapCalculation(imagePyr_3, imagePyr_3);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_4, featureType);
	BlockDivision(imagePyr_4, imagePyr_4);
	GlobalMapCalculation(imagePyr_4, imagePyr_4);

	PyramidSummation(dst);
}

void GlobalCues::GlobalMapCalculation(InputArray src, OutputArray dst) {
	Mat srcMat = src.getMat();
	int x = srcMat.rows, y = srcMat.cols;

	if(srcMat.channels() > 1) {
		Mat resultMat(x, y, CV_32FC3);

		for(int i = 0; i < x; i++){
			for(int j = 0; j < y; j++) {
				Mat diffMat(x, y, srcMat.type());
				//int dep = srcMat.depth();
				double pixelValue_0 = srcMat.at<Vec3b>(i,j)[0];
				double pixelValue_1 = srcMat.at<Vec3b>(i,j)[1];
				double pixelValue_2 = srcMat.at<Vec3b>(i,j)[2];
				Scalar pixelValue;
				pixelValue.val[0] = pixelValue_0;
				pixelValue.val[1] = pixelValue_1;
				pixelValue.val[2] = pixelValue_2;
				pixelValue.val[3] = 0.0;
				absdiff(srcMat, pixelValue, diffMat);
				Scalar meanValue = mean(diffMat);
				resultMat.at<Vec3f>(i,j)[0] = meanValue.val[0];
				resultMat.at<Vec3f>(i,j)[1] = meanValue.val[1];
				resultMat.at<Vec3f>(i,j)[2] = meanValue.val[2];
				diffMat.release();
			}
		}
		resultMat.copyTo(dst);
		resultMat.release();
	}
	else {
		Mat resultMat(x, y, DataType<uchar>::type);

		for(int i = 0; i < x; i++){
			for(int j = 0; j < y; j++) {
				Mat diffMat(x, y, srcMat.type());
				uchar pixelValue = srcMat.at<uchar>(i,j);
				absdiff(srcMat, Scalar::all(pixelValue), diffMat);
				Scalar meanValue = mean(diffMat);
				resultMat.at<uchar>(i,j) = meanValue.val[0];
				diffMat.release();
			}
		}

		normalize(resultMat, resultMat, 255, 0, NORM_L2, -1, noArray());
		resultMat.copyTo(dst);
		resultMat.release();
	}
}
