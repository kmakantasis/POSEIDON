/*
 * Classifier.cpp
 *
 *  Created on: Feb 8, 2013
 *      Author: kostas makantasis
 */

#include "Classifier.h"

Classifier::Classifier() {}

Classifier::~Classifier() {}

Mat Classifier::CreateFeatureVec(InputArray src) {

	Mat srcMat = src.getMat();
	cvtColor(srcMat, srcMat, CV_RGB2GRAY);
	Mat dstMat;

	dstMat = srcMat.reshape(0, srcMat.rows*srcMat.cols);

	return dstMat;
}

void Classifier::BayesTrain(const Mat& src, const Mat& responses) {

	classifier.train(src, responses, Mat(), Mat(), false);
}

void Classifier::BayesPredict(const Mat& src, Mat* results) {

	classifier.predict(src, results);
}
