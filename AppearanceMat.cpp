/*
 * AppearanceMat.cpp
 *
 *  Created on: Feb 8, 2013
 *      Author: kostas makantasis
 */

#include "AppearanceMat.h"

AppearanceMat::AppearanceMat() {}

AppearanceMat::~AppearanceMat() {}

void AppearanceMat::SingleCueMat(InputArray src1, InputArray src2, InputArray src3, InputArray src4, OutputArray dst) {

	vector<Mat> colorPlanes;

	addWeighted(src1, 0.25, src2, 0.25, 0, dst);
	addWeighted(dst, 1.0, src3, 0.25, 0, dst);
	split(src4, colorPlanes);
	addWeighted(dst, 1.0, colorPlanes[0], 0.08, 0, dst);
	addWeighted(dst, 1.0, colorPlanes[1], 0.08, 0, dst);
	addWeighted(dst, 1.0, colorPlanes[2], 0.08, 0, dst);

	//applyColorMap(dst, dst, COLORMAP_JET);
}

void AppearanceMat::AllCuesMat(InputArray src1, InputArray src2, InputArray src3, OutputArray dst) {

	addWeighted(src1, 0.33, src2, 0.33, 0, dst);
	addWeighted(dst, 1.0, src3, 0.33, 0, dst);
	applyColorMap(dst, dst, COLORMAP_JET);

}
