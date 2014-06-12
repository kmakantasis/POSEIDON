/*
 * CSCues.cpp
 *
 *  Created on: Feb 2, 2013
 *      Author: kostas makantasis
 */

#include "CSCues.h"

CSCues::CSCues(double div) {this->blockDiv = div;}

CSCues::~CSCues() {}

void CSCues::CSCuesImage(InputArray src, OutputArray dst, int featureType) {

	Mat srcMat = src.getMat();
	original = srcMat.clone();

	PyramidFeature(src, imagePyr_1, featureType);
	BlockDivision(imagePyr_1, imagePyr_1);
	CSMapCalculation(imagePyr_1, imagePyr_1);

	CreatePyramid(src, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_2, featureType);
	BlockDivision(imagePyr_2, imagePyr_2);
	CSMapCalculation(imagePyr_2, imagePyr_2);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_3, featureType);
	BlockDivision(imagePyr_3, imagePyr_3);
	CSMapCalculation(imagePyr_3, imagePyr_3);

	CreatePyramid(pyrTemp, pyrTemp);
	PyramidFeature(pyrTemp, imagePyr_4, featureType);
	BlockDivision(imagePyr_4, imagePyr_4);
	CSMapCalculation(imagePyr_4, imagePyr_4);

	PyramidSummation(dst);
}

void CSCues::CSMapCalculation(InputArray src, OutputArray dst){

	Mat srcMat = src.getMat();
	int windowSizeX, windowSizeY;
	int x = srcMat.rows, y = srcMat.cols;
	Point point1, point2;
	Rect rect;

	if(srcMat.channels() > 1) {
			Mat resultMat(x, y, CV_32FC3);
			for(int i = 0; i < x; i++){
				if(i < x-i)
					windowSizeX = i;
				else
					windowSizeX = x-i;
				for(int j = 0; j < y; j++) {
					if(j < y-j)
						windowSizeY = j;
					else
						windowSizeY = y-j;
					point1.y = i-windowSizeX, point1.x = j-windowSizeY;
					point2.y = i+windowSizeX, point2.x = j+windowSizeY;
					rect = Rect(point1.x,point1.y,point2.x-point1.x,point2.y-point1.y);
					Mat roiMat = srcMat(rect);
					Mat diffMat(roiMat.rows, roiMat.cols, roiMat.type());
					double pixelValue_0 = srcMat.at<Vec3b>(i,j)[0];
					double pixelValue_1 = srcMat.at<Vec3b>(i,j)[1];
					double pixelValue_2 = srcMat.at<Vec3b>(i,j)[2];
					Scalar pixelValue;
					pixelValue.val[0] = pixelValue_0;
					pixelValue.val[1] = pixelValue_1;
					pixelValue.val[2] = pixelValue_2;
					pixelValue.val[3] = 0.0;
					absdiff(roiMat, pixelValue, diffMat);
					Scalar meanValue = mean(diffMat);
					resultMat.at<Vec3f>(i,j)[0] = meanValue.val[0];
					resultMat.at<Vec3f>(i,j)[1] = meanValue.val[1];
					resultMat.at<Vec3f>(i,j)[2] = meanValue.val[2];
					diffMat.release();
					roiMat.release();
				}
			}
			resultMat.copyTo(dst);
			resultMat.release();
	}
	else {
		Mat resultMat(x, y, DataType<uchar>::type);
		for(int i = 0; i < x; i++){
			if(i < x-i)
				windowSizeX = i;
			else
				windowSizeX = x-i;
			for(int j = 0; j < y; j++) {
				if(j < y-j)
					windowSizeY = j;
				else
					windowSizeY = y-j;
				point1.y = i-windowSizeX, point1.x = j-windowSizeY;
				point2.y = i+windowSizeX, point2.x = j+windowSizeY;
				rect = Rect(point1.x,point1.y,point2.x-point1.x,point2.y-point1.y);
				Mat roiMat = srcMat(rect);
				Mat diffMat(roiMat.rows, roiMat.cols, roiMat.type());
				uchar pixelValue = srcMat.at<uchar>(i,j);
				absdiff(roiMat, Scalar::all(pixelValue), diffMat);
				Scalar meanValue = mean(diffMat);
				resultMat.at<uchar>(i,j) = meanValue.val[0];
				diffMat.release();
				roiMat.release();
			}
		}
		normalize(resultMat, resultMat, 255, 0, NORM_L2, -1, noArray());
		resultMat.copyTo(dst);
		resultMat.release();
	}

}
