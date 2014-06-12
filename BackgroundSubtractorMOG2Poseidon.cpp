/*
 * BackgroundSubtractorMOG2Poseidon.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: kostas
 */

#include "BackgroundSubtractorMOG2Poseidon.h"

BackgroundSubtractorMOG2_Poseidon::BackgroundSubtractorMOG2_Poseidon() {}

BackgroundSubtractorMOG2_Poseidon::BackgroundSubtractorMOG2_Poseidon(int history, float varThreshold, bool bShadowDetection) {
			BackgroundSubtractorMOG2(history, varThreshold, bShadowDetection);
}

BackgroundSubtractorMOG2_Poseidon::~BackgroundSubtractorMOG2_Poseidon() {}

void BackgroundSubtractorMOG2_Poseidon::setHistory(int k) {
	history = k;
}

void BackgroundSubtractorMOG2_Poseidon::setNmixtures(int k) {
	nmixtures = k;
}

void BackgroundSubtractorMOG2_Poseidon::setVarThreshold(double k) {
	varThreshold = k;
}

void BackgroundSubtractorMOG2_Poseidon::setBackgroundRatio(float k) {
	backgroundRatio = k;
}

void BackgroundSubtractorMOG2_Poseidon::setVarThresholdGen(float k) {
	varThresholdGen = k;
}

void BackgroundSubtractorMOG2_Poseidon::setFVarInit(float k) {
	fVarInit = k;
}

void BackgroundSubtractorMOG2_Poseidon::setFVarMin(float k) {
	fVarMin = k;
}

void BackgroundSubtractorMOG2_Poseidon::setFVarMax(float k) {
	fVarMax = k;
}

void BackgroundSubtractorMOG2_Poseidon::setFCT(float k) {
	fCT = k;
}

void BackgroundSubtractorMOG2_Poseidon::setBShadowDetection(bool k) {
	bShadowDetection = k;
}

void BackgroundSubtractorMOG2_Poseidon::setNShadowDetection(uchar k) {
	nShadowDetection = k;
}

void BackgroundSubtractorMOG2_Poseidon::setFTau(float k) {
	fTau = k;
}

