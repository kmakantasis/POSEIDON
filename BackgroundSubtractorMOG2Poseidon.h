/*
 * BackgroundSubtractorMOG2Poseidon.h
 *
 *  Created on: Feb 11, 2013
 *      Author: kostas
 */

#ifndef BACKGROUNDSUBTRACTORMOG2POSEIDON_H_
#define BACKGROUNDSUBTRACTORMOG2POSEIDON_H_

#include <opencv2/opencv.hpp>

class BackgroundSubtractorMOG2_Poseidon: public cv::BackgroundSubtractorMOG2 {
public:
	BackgroundSubtractorMOG2_Poseidon();
	BackgroundSubtractorMOG2_Poseidon(int history, float varThreshold, bool bShadowDetection);
	virtual ~BackgroundSubtractorMOG2_Poseidon();
	void setHistory(int k);
	void setNmixtures(int k);
	void setVarThreshold(double k);
	void setBackgroundRatio(float k);
	void setVarThresholdGen(float k);
	void setFVarInit(float k);
	void setFVarMin(float k);
	void setFVarMax(float k);
	void setFCT(float k);
	void setBShadowDetection(bool k);
	void setNShadowDetection(uchar k);
	void setFTau(float k);
};

#endif /* BACKGROUNDSUBTRACTORMOG2POSEIDON_H_ */
