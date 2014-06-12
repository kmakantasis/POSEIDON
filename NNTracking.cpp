/*
 * NNTracking.cpp
 *
 *  Created on: Mar 9, 2013
 *      Author: kostas makantasis
 */

#include "NNTracking.h"

NNTracking::NNTracking() {}

NNTracking::~NNTracking() {}

void NNTracking::CreateNN() {

	inputSamples = imread("../Data/inputSamples.jpg");
	cvtColor(inputSamples, inputSamples, CV_RGB2GRAY);
	inputSamples.convertTo(inputSamples, CV_32FC1);
	outputSamples = imread("../Data/outputSamples2.jpg");
	cvtColor(outputSamples, outputSamples, CV_RGB2GRAY);
	outputSamples.convertTo(outputSamples, CV_32FC1);

	Mat layers = Mat(1,3,CV_32SC1);
	layers.at<int>(0,0) = inputSamples.cols;
	layers.at<int>(0,1) = 10;
	layers.at<int>(0,2) = 2	;

	nnetwork.create(layers, CvANN_MLP::SIGMOID_SYM, 0.6, 1);

}

void NNTracking::TrainNN() {
	 CvANN_MLP_TrainParams params = CvANN_MLP_TrainParams(cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 1000, 0.000001), CvANN_MLP_TrainParams::BACKPROP, 0.1, 0.1);

	 int iterations = nnetwork.train(inputSamples, outputSamples, Mat(), Mat(), params);

	 cout << "Training iterations: " << iterations << endl;
}

Mat NNTracking::PredictNN(InputArray inputs) {

	Mat nnIn = inputs.getMat();
	int t = nnIn.depth();
	Mat nnOut;

	nnetwork.predict(nnIn, nnOut);

	nnIn.release();

	return nnOut;

}

void NNTracking::Retrain(InputArray newTrainingInput, InputArray newTrainingOutput) {
	Mat newInputs;
	Mat newOutputs;

	newInputs = newTrainingInput.getMat();
	newOutputs = newTrainingOutput.getMat();

	Mat layers = Mat(1,3,CV_32SC1);
	layers.at<int>(0,0) = inputSamples.cols;
	layers.at<int>(0,1) = 10;
	layers.at<int>(0,2) = 2	;

	nnetworkRetrained.create(layers, CvANN_MLP::SIGMOID_SYM, 0.6, 1);
	CvANN_MLP_TrainParams params = CvANN_MLP_TrainParams(cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 1000, 0.000001), CvANN_MLP_TrainParams::BACKPROP, 0.1, 0.1);
	int iterations = nnetwork.train(newInputs, newOutputs, Mat(), Mat(), params, CvANN_MLP::UPDATE_WEIGHTS);

}











