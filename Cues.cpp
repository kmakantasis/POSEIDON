/*
 * Cues.cpp
 *
 *  Created on: Jan 30, 2013
 *      Author: kostas makantasis
 */

#include "Cues.h"
#include "FrequencyFeatures.h"
#include "EdgeFeatures.h"
#include "RightAngleFeatures.h"
#include "ColorFeatures.h"

Cues::Cues() {}

Cues::~Cues() {}

void Cues::BlockDivision(InputArray src, OutputArray dst){

	resize(src, dst, Size(), blockDiv, blockDiv, INTER_AREA);

}

void Cues::CreatePyramid(InputArray src, OutputArray dst){

	pyrDown(src, dst, Size());

}

int Cues::PyramidFeature(InputArray src, OutputArray dst, int featureType){

	if(featureType == 1) { // Edge feature selected
		EdgeFeatures edges;
		edges.FindEdges(src, dst, 25, 3, 3);
	}
	if(featureType == 2) { // Frequency feature selected
		FrequencyFeatures frequencies;
		frequencies.FindFrequency(src, dst);
	}
	if(featureType == 3) { // Right angles feature selected
		RightAngleFeatures rightAngles;
		rightAngles.FindRightAngles(src, dst);
	}
	if(featureType == 41 || featureType == 42) { // Color feature selected
		ColorFeatures colorFeat;
		colorFeat.FindColor(src, dst, featureType);
	}


	return 0;
}

void Cues::ColorMap(InputArray src, OutputArray dst) {

	resize(src, dst, original.size(), 0, 0);
	//applyColorMap(dst, dst, COLORMAP_JET);
	//normalize(dst, dst, 255, 0, NORM_L2, CV_32F, noArray());
}

void Cues::PyramidSummation(OutputArray dst){

	resize(imagePyr_4, scaleTemp, imagePyr_3.size(), 0, 0);
	add(imagePyr_3, scaleTemp, sumTempPyr34, noArray(), imagePyr_3.depth());

	resize(imagePyr_3, scaleTemp, imagePyr_2.size(), 0, 0);
	add(imagePyr_2, scaleTemp, sumTempPyr23);

	resize(imagePyr_2, scaleTemp, imagePyr_1.size(), 0, 0);
	add(imagePyr_1, scaleTemp, sumTempPyr12);

	resize(sumTempPyr34, scaleTemp, sumTempPyr23.size(), 0, 0);
	add(sumTempPyr23, scaleTemp, sumTemp);

	resize(sumTemp, scaleTemp, sumTempPyr12.size(), 0, 0);
	add(sumTempPyr12, scaleTemp, sumTemp);

	double alpha = 1;

	convertScaleAbs(sumTemp, dst, alpha);

	ColorMap(dst, dst);

}
