/*
 * AppearanceMat.h
 *
 *  Created on: Feb 8, 2013
 *      Author: kostas makantasis
 */

#include <opencv2/opencv.hpp>

using namespace cv;


#ifndef APPEARANCEMAT_H_
#define APPEARANCEMAT_H_

class AppearanceMat {
public:
	AppearanceMat();
	virtual ~AppearanceMat();
	void SingleCueMat(InputArray src1, InputArray src2, InputArray src3, InputArray src4, OutputArray dst);
	void AllCuesMat(InputArray src1, InputArray src2, InputArray src3, OutputArray dst);
};

#endif /* APPEARANCEMAT_H_ */
