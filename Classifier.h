/*
 * Classifier.h
 *
 *  Created on: Feb 8, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;

#ifndef CLASSIFIER_H_
#define CLASSIFIER_H_

class Classifier {
protected:
	CvNormalBayesClassifier classifier;
public:
	Classifier();
	virtual ~Classifier();
	Mat CreateFeatureVec(InputArray src);
	void BayesTrain(const Mat& src, const Mat& responses);
	void BayesPredict(const Mat& src, Mat* results);
};

#endif /* CLASSIFIER_H_ */
