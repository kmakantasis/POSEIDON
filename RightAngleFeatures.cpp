/*
 * RightAngleFeatures.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: kostas makantasis
 */

#include "RightAngleFeatures.h"

RightAngleFeatures::RightAngleFeatures() {}

RightAngleFeatures::~RightAngleFeatures() {}

void RightAngleFeatures::FindRightAngles(InputArray src, OutputArray dst) {

	EdgeFeatures edges;
	edges.FindEdges(src, imageEdges, 25, 3, 3);


	/***** KERNEL INITIALIZATION ********
	 *      |0   0   1   0   0|         *
	 *      |0   0   2   0   0|			*
	 *  k = |1   2   4   2   1| * 1/16	*
	 *      |0   0   2   0   0|			*
	 *      |0   0   1   0   0|			*
	 ************************************/

	Mat kernel;
	int kernel_size = 5;
	kernel = Mat::zeros( kernel_size, kernel_size, CV_32F );

	kernel.at<float>(0,2) = 1 / (float)(kernel_size*kernel_size);
	kernel.at<float>(1,2) = 2 / (float)(kernel_size*kernel_size);
	kernel.at<float>(2,0) = 1 / (float)(kernel_size*kernel_size);
	kernel.at<float>(2,1) = 2 / (float)(kernel_size*kernel_size);
	kernel.at<float>(2,2) = 4 / (float)(kernel_size*kernel_size);
	kernel.at<float>(2,3) = 2 / (float)(kernel_size*kernel_size);
	kernel.at<float>(2,4) = 1 / (float)(kernel_size*kernel_size);
	kernel.at<float>(3,2) = 2 / (float)(kernel_size*kernel_size);
	kernel.at<float>(4,2) = 1 / (float)(kernel_size*kernel_size);

	filter2D(imageEdges, dst, -1 , kernel, Point(-1, -1), 0, BORDER_DEFAULT );

}
