/*
 * NNTracking.h
 *
 *  Created on: Mar 9, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>
#include "CvANNMLPExt.h"

using namespace cv;
using namespace std;

#ifndef NNTRACKING_H_
#define NNTRACKING_H_

class NNTracking {
protected:
	CvANN_MLP_Ext nnetwork;
	CvANN_MLP_Ext nnetworkRetrained;
	Mat inputSamples;
	Mat outputSamples;
	double* weightsL0;
	double* weightsL1;
	double* weightsL2;
	double* retWeightsL0;
	double* retWeightsL1;
	double* retWeightsL2;
	double* newWeight;
public:
	NNTracking();
	virtual ~NNTracking();
	void CreateNN();
	void TrainNN();
	void Retrain(InputArray newTrainingInput, InputArray newTrainingOutput);
	Mat PredictNN(InputArray inputs);
};

#endif /* NNTRACKING_H_ */
