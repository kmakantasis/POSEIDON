/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: kostas makantasis
 */
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <string.h>
//#include <curl/curl.h>
#include "EdgeFeatures.h"
#include "LocalCues.h"
#include "GlobalCues.h"
#include "CSCues.h"
#include "AppearanceMat.h"
#include "BackgroundSubtractorMOG2Poseidon.h"
#include "poseidontest.h"
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <string>
#include "config.h"
#include "foregroundinterface.h"
#include "videocapture.h"
#include "foreground.h"
#include "datasetinfo.h"
#include "NNTracking.h"

#define EDGES 1
#define FREQUENCY 2
#define RIGHTANGLES 3
#define COLOR_LOCAL 41
#define COLOR_GCS 42

using namespace cv;
using namespace std;

//////// TRACKBAR VARIABLES /////////
const int alpha_slider_max = 100;
int alpha_slider;
double alpha;
double beta;
Mat trackSrc1, trackSrc2, trackDst;
/////////////////////////////////////

//////// TRACKBAR VARIABLES /////////
const int wave_slider_max = 100;
int wave_slider;
double wave;
/////////////////////////////////////


void on_trackbar(int, void*) {
	alpha = (double) alpha_slider/alpha_slider_max ;
	beta = ( 1.0 - alpha );

	addWeighted( trackSrc1, alpha, trackSrc2, beta, 0.0, trackDst);
}

void on_trackbar_waves(int, void*) {
	wave = (double) wave_slider + 100;
}

Mat BGSubtraction(BackgroundSubtractorMOG2_Poseidon* bg, InputArray src, InputArray prob) {
	Mat frame, frame_bw, frame_blur, back, back_bw, fore, back_mask, back_mask2;
	std::vector<std::vector<cv::Point> > contours;
	vector<Mat> redPlane;

	frame = src.getMat();
	split(prob, redPlane);

    GaussianBlur(frame,frame_blur,cvSize(5,5),2.5);
    cvtColor(frame_blur, frame_bw, CV_RGB2GRAY);
    bg->operator()(frame_blur,fore);  //fore-->foregroundmask
    bg->getBackgroundImage(back);
    cvtColor(back, back_bw, CV_RGB2GRAY);
    subtract(back_bw,frame_bw,back_mask);
    threshold(back_mask, back_mask, 25, 255, THRESH_BINARY);
    convertScaleAbs(redPlane[0], redPlane[0], 1, -255);
    //threshold(redPlane[0], redPlane[0], 128, 255, THRESH_BINARY_INV);
    dilate(redPlane[0],redPlane[0],cv::Mat());
    dilate(redPlane[0],redPlane[0],cv::Mat());

    medianBlur(fore,fore,7);
    medianBlur(fore,fore,7);
    erode(fore,fore,cv::Mat());
    erode(fore,fore,cv::Mat());
    erode(fore,fore,cv::Mat());
    dilate(fore,fore,cv::Mat());
    dilate(fore,fore,cv::Mat());
    dilate(fore,fore,cv::Mat());
    dilate(fore,fore,cv::Mat());
    //medianBlur(back_mask,back_mask,3);
    dilate(back_mask,back_mask,cv::Mat());
    dilate(back_mask,back_mask,cv::Mat());
    absdiff(fore,back_mask,back_mask2);
    threshold(back_mask2, back_mask2, 25, 255, THRESH_BINARY);
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());
    dilate(back_mask2,back_mask2,cv::Mat());

    //////////// TRACKBAR /////////////
    trackSrc1 = back_mask2.clone();
    trackSrc2 = redPlane[0].clone();
    on_trackbar( alpha_slider,0);
    back_mask2 = trackDst.clone();

    on_trackbar_waves( wave_slider,0);
    ///////////////////////////////////

    threshold(back_mask2, back_mask2, wave, 255, THRESH_BINARY);

    findContours(back_mask2,contours,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE);
    drawContours(frame,contours,-1,cv::Scalar(0,0,255),1);

    return frame;
}


Mat NNOutput(InputArray salientOutput, InputArray bgOutput, NNTracking* NN, int r) {
	Mat nnOutput;
	Mat fnnOutput;
	Mat nnOutput2C;
	Mat newInput;
	Mat bgTemp;
	Mat sTemp;
	Mat newTrainingInput(750,13,CV_32FC1);
	Mat newTrainingOutput(750,2,CV_32FC1);
	Mat difference;

	sTemp = salientOutput.getMat();
	bgTemp = bgOutput.getMat();
	bgTemp = bgTemp.reshape(0, bgTemp.rows*bgTemp.cols);
	//cvtColor(bgTemp, bgTemp, CV_RGB2GRAY);

	int sType = sTemp.type();
	int bgType = bgTemp.type();

	hconcat(sTemp, bgTemp, newInput);

	newInput.convertTo(newInput, CV_32FC1);

	int newType = newInput.type();

	nnOutput2C = NN->PredictNN(newInput);
	nnOutput = nnOutput2C.col(0).clone();
	fnnOutput = nnOutput.reshape(0, r);

	//Mat aDiff;
	//aDiff.convertTo(fnnOutput, CV_8U, 255);

	//absdiff(aDiff, bgOutput, difference);

	int count = 0, flagF = 0, flagB = 0;
	for(int i=0; i<nnOutput2C.rows;i++){
		if(nnOutput.at<float>(i,0) > 0.9){
			newInput.row(i).copyTo(newTrainingInput.row(count));
			newTrainingOutput.at<float>(count,0) = 0.0;
			newTrainingOutput.at<float>(count,1) = 1.0;
			count++;
		}
		if(count == 500) {
			flagB = 1;
			break;
		}
	}
	for(int i=0; i<nnOutput2C.rows;i++){
		if(nnOutput.at<float>(i,0) < 0.1){
			newInput.row(i).copyTo(newTrainingInput.row(count));
			newTrainingOutput.at<float>(count,0) = 1.0;
			newTrainingOutput.at<float>(count,1) = 0.0;
			count++;
		}
		if(count == 750) {
			flagF = 1;
			break;
		}
	}
	/////////////// New training set for retraining is ready ////////////////
	if (flagF ==1 && flagB ==1)
		NN->Retrain(newTrainingInput, newTrainingOutput);
	///////////////  I have to call the Retrain() function   ////////////////
	else {
		threshold(fnnOutput, fnnOutput, 0.95, 1, THRESH_BINARY);
		imshow( "Network Output", fnnOutput);
		imwrite( "results4/fnnOutput.jpg", fnnOutput);
	}
	cout << fnnOutput.at<float>(39,273);
	cout << endl;
	return fnnOutput;
}


Mat ProbabilityMap(InputArray src, LocalCues l, GlobalCues g, CSCues c, AppearanceMat a, OutputArray salientOutput) {

	Mat image, resultLocal, resultGlobal, resultCS, resultLocalF, resultGlobalF, resultCSF, resultLocalR, resultGlobalR, resultCSR;
	Mat resultLocalC, resultGlobalC, resultCSC, resultLocalE, resultGlobalE, resultCSE, resultTotal;
	image = src.getMat();

	l.LocalCuesImage(image, resultLocalE, EDGES);
	g.GlobalCuesImage(image, resultGlobalE, EDGES);
	c.CSCuesImage(image, resultCSE, EDGES);

	l.LocalCuesImage(image, resultLocalF, FREQUENCY);
	g.GlobalCuesImage(image, resultGlobalF, FREQUENCY);
	c.CSCuesImage(image, resultCSF, FREQUENCY);

	l.LocalCuesImage(image, resultLocalR, RIGHTANGLES);
	g.GlobalCuesImage(image, resultGlobalR, RIGHTANGLES);
	c.CSCuesImage(image, resultCSR, RIGHTANGLES);

	l.LocalCuesImage(image, resultLocalC, COLOR_LOCAL);
	g.GlobalCuesImage(image, resultGlobalC, COLOR_GCS);
	c.CSCuesImage(image, resultCSC, COLOR_GCS);

	a.SingleCueMat(resultLocalE, resultLocalF, resultLocalR, resultLocalC, resultLocal);
	a.SingleCueMat(resultGlobalE, resultGlobalF, resultGlobalR, resultGlobalC, resultGlobal);
	a.SingleCueMat(resultCSE, resultCSF, resultCSR, resultCSC, resultCS);

	/*imwrite( "edges_local.jpg", resultLocalE );
	imwrite( "edges_global.jpg", resultGlobalE );
	imwrite( "edged_neigh.jpg", resultCSE );
	imwrite( "frequency_local.jpg", resultLocalF );
	imwrite( "frequency_global.jpg", resultGlobalF );
	imwrite( "frequency_neigh.jpg", resultCSF );
	imwrite( "lines_local.jpg", resultLocalR );
	imwrite( "lines_global.jpg", resultGlobalR );
	imwrite( "lines_neigh.jpg", resultCSR );
	imwrite( "color_local.jpg", resultLocalC );
	imwrite( "color_global.jpg", resultGlobalC );
	imwrite( "color_neigh.jpg", resultCSC );*/

	a.AllCuesMat(resultLocal, resultGlobal, resultCS, resultTotal);

	////////////////////////////////////////////////////////////////////
	Mat le = resultLocalE.reshape(0,resultLocalE.rows*resultLocalE.cols);
	Mat ge = resultGlobalE.reshape(0,resultGlobalE.rows*resultGlobalE.cols);
	Mat cse = resultCSE.reshape(0,resultCSE.rows*resultCSE.cols);

	Mat lf = resultLocalE.reshape(0,resultLocalF.rows*resultLocalF.cols);
	Mat gf = resultGlobalE.reshape(0,resultGlobalF.rows*resultGlobalF.cols);
	Mat csf = resultCSE.reshape(0,resultCSF.rows*resultCSF.cols);

	Mat lr = resultLocalR.reshape(0,resultLocalR.rows*resultLocalR.cols);
	Mat gr = resultGlobalR.reshape(0,resultGlobalR.rows*resultGlobalR.cols);
	Mat csr = resultCSR.reshape(0,resultCSR.rows*resultCSR.cols);

	cvtColor(resultLocalC, resultLocalC, CV_RGB2GRAY);
	cvtColor(resultGlobalC, resultGlobalC, CV_RGB2GRAY);
	cvtColor(resultCSC, resultCSC, CV_RGB2GRAY);

	Mat lc = resultLocalC.reshape(0,resultLocalC.rows*resultLocalC.cols);
	Mat gc = resultGlobalC.reshape(0,resultGlobalC.rows*resultGlobalC.cols);
	Mat csc = resultCSC.reshape(0,resultCSC.rows*resultCSC.cols);

	hconcat(lr, gr, salientOutput);
	hconcat(salientOutput, csr, salientOutput);
	hconcat(salientOutput, lf, salientOutput);
	hconcat(salientOutput, gf, salientOutput);
	hconcat(salientOutput, csf, salientOutput);
	hconcat(salientOutput, le, salientOutput);
	hconcat(salientOutput, ge, salientOutput);
	hconcat(salientOutput, cse, salientOutput);
	hconcat(salientOutput, lc, salientOutput);
	hconcat(salientOutput, gc, salientOutput);
	hconcat(salientOutput, csc, salientOutput);

	return resultTotal;

}

/*void PTZControl(char* action) {
	CURL *curl;
	CURLcode res;

	curl = curl_easy_init();
	if(curl) {
		curl_easy_setopt(curl, CURLOPT_URL, action);
	    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

	    res = curl_easy_perform(curl);
	    if(res != CURLE_OK)
	    	fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));

	    curl_easy_cleanup(curl);
	  }
}*/

Mat IterativeFunc(Mat image, LocalCues localCues, GlobalCues globalCues, CSCues csCues, AppearanceMat ApMat, Mat totalPost, BackgroundSubtractorMOG2_Poseidon* bg, OutputArray dst, OutputArray salientOutput, OutputArray netOutput, NNTracking* NN, int bgFlag) {
	clock_t tStart = clock();
	Mat original = image.clone();
	Mat prob = ProbabilityMap(image, localCues, globalCues, csCues, ApMat, salientOutput);
	Mat frame;
	if (bgFlag == 0){
		totalPost = BGSubtraction(bg, image, prob);
		cvtColor(totalPost, totalPost, CV_RGB2GRAY);
		//totalPost = frame.clone() //////////////////////////////////////////////////////
	}
	totalPost.convertTo(totalPost, CV_8U, 255);
	int totalPostType = totalPost.type();
	int totalPostChannels = totalPost.channels();

	Mat nnOutput = NNOutput(salientOutput, totalPost, NN, totalPost.rows);

	imshow( "Display Original", original );
	imshow( "Display BG Mask", totalPost );
	imshow( "Display Probability", prob );


	cvWaitKey(1);

	prob.copyTo(dst);
	nnOutput.copyTo(netOutput);

	return original;

}

#define INTEGRATION
#ifdef INTEGRATION
t_settings* g_pSetts;
t_horizon* g_pHorizon;

//int mainPoseidon(int argc, char *argv[])
int mainPoseidon(char *input_file, char *output_file)
{
	t_info info, *pInfo=&info;
	t_vid vid, *pVid=&vid;
	t_bg bg, *pBg=&bg;
	t_horizon horizon, *pHorizon=&horizon;
	g_pHorizon = pHorizon;
	t_track aTrack[MAX_TARGETS];
	t_settings settings;
	g_pSetts=&settings;
	std::string szPath;

	//getPath(argv[1], szPath);
	getPath(input_file, szPath);
	readConfigFile(szPath, g_pSetts);
	createDataSetInfo(pInfo);
	createHorizon(g_pSetts->U, g_pSetts->V, pHorizon);

	//std::string szFilename = argv[1];
	std::string szFilename = input_file;
	std::string szFilenameWithPath = szPath + szFilename ;
	createVideoCapture(szPath, szFilenameWithPath, g_pSetts->U, g_pSetts->V, pVid);
	initialiseEverythingParis(pInfo, pVid, pHorizon, pBg, aTrack);


	//int k = PoseidonTest(argc, argv[1]);

	Mat image, frame, prob, original, nnOutput;
	Mat salientOutput;
	double totalTime;
	int framesSync, wKey, totalTimeInt;
	string filename;
	//char c;
	//char* action;


	///************************MOG INITIALIATION*****************************///
	BackgroundSubtractorMOG2_Poseidon* kmBg = new BackgroundSubtractorMOG2_Poseidon();
	kmBg->setHistory(50);
	kmBg->setNmixtures(5);
	kmBg->setVarThreshold(9.0); kmBg->setVarThresholdGen(9.0);
	kmBg->setBackgroundRatio(0.6);
	kmBg->setFVarInit(15.0);
	kmBg->setFCT(0.01);
	kmBg->setBShadowDetection(true); kmBg->setNShadowDetection(127); kmBg->setFTau(0.5);
	////////////////////////////////////////////////////////////////////////////

	NNTracking NN;
	NN.CreateNN();
	NN.TrainNN();



	//clock_t tStart = clock();
	//VideoCapture cap("http://root:poseidon@147.27.11.212:80/mjpg/video.mjpg");
	//if(strcmp(argv[1],"web") == 0) {
	if(strcmp(input_file,"web") == 0) {
		filename = "http://root:poseidon@147.27.11.212:80/mjpg/video.mjpg";
	}
	else {
		//filename = argv[1];
		filename = input_file;
	}


	VideoCapture cap(filename);
	cap >> image;

	double blockSize = (double)80/image.cols;
	LocalCues localCues(blockSize);
	GlobalCues globalCues(blockSize);
	CSCues csCues(blockSize);
	AppearanceMat ApMat;

	namedWindow( "Display Original", CV_WINDOW_AUTOSIZE );
	namedWindow( "Display BG Mask", CV_WINDOW_AUTOSIZE );
	namedWindow( "Display Probability", CV_WINDOW_AUTOSIZE );


	//////////////////////////// TRACKBAR INITIALIZATION ////////////////////////////////////////////
	alpha_slider = 40;
	char TrackbarName[50];
	sprintf( TrackbarName, "BG Infl.");
	createTrackbar( TrackbarName, "Display BG Mask", &alpha_slider, alpha_slider_max, on_trackbar );

	char WaveTrackbarName[50];
	sprintf( WaveTrackbarName, "Waves");
	wave_slider = 25;
	createTrackbar( WaveTrackbarName, "Display BG Mask", &wave_slider, wave_slider_max, on_trackbar_waves );
	/////////////////////////////////////////////////////////////////////////////////////////////////

	Mat prob32bit, totalPostTemp;
	bool stop = 0;

	//VideoWriter output_cap(output_file,
	//               cap.get(CV_CAP_PROP_FOURCC),
	//               cap.get(CV_CAP_PROP_FPS),
	//               cv::Size(cap.get(CV_CAP_PROP_FRAME_WIDTH),
	//               cap.get(CV_CAP_PROP_FRAME_HEIGHT)));

	int bgFlag = 0;
	Mat sumPost, totalPost;
	int count = 0;

	while(!stop) {
		cap >> image;
		std::ostringstream save_orig;
		save_orig << "results95/original" << count << ".jpg";
		imwrite( save_orig.str(), image);
		Mat orig;
		orig = image.clone();
		IplImage converted = IterativeFunc(image, localCues, globalCues, csCues, ApMat, totalPost, kmBg, prob, salientOutput, nnOutput, &NN, bgFlag);



		cvResize(&converted, pVid->Iresized);

		int bContinue = doEverything(pInfo, pVid, pHorizon, pBg, aTrack);
		Mat posteriorBG(pBg->IjointBG);
		resize(posteriorBG, posteriorBG, Size(image.cols,image.rows));
		imshow("posteriorBG", posteriorBG);
		std::ostringstream save_posterioBG;
		save_posterioBG << "results95/posteriorBG" << count << ".jpg";
		imwrite( save_posterioBG.str(), posteriorBG);


		imshow("posteriorKM", prob);
		std::ostringstream save_posterioKM;
		save_posterioKM << "results95/posteriorKM" << count << ".jpg";
		imwrite( save_posterioKM.str(), prob);
		cvtColor(prob, prob, CV_BGR2GRAY);

		prob.convertTo(prob32bit, CV_32FC1, 1/255.0);
		double minV, maxV;
		minMaxLoc(prob32bit, &minV, &maxV);



		cout << prob32bit.cols <<"x"<< prob32bit.rows << " " << prob32bit.channels() << endl;
		cout << posteriorBG.cols <<"x"<< posteriorBG.rows << endl;
		cout.flush();

		add(prob32bit, posteriorBG, sumPost);
		divide(prob32bit, sumPost, totalPost);

		Mat postThreshold;
		threshold(totalPost, postThreshold, 0.95, 255, THRESH_BINARY);

		//totalPost.convertTo(totalPost, CV_8UC1, 255);
		//std::vector<std::vector<cv::Point> > contours;
		//cout << totalPost.type() << " " << totalPost.channels() << endl;
		//findContours(totalPost,contours,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE);
		//drawContours(orig,contours,-1,cv::Scalar(0,0,255),2);

		imshow("Total Posterior", postThreshold);
		std::ostringstream save_totalPosterior;
		save_totalPosterior << "results95/totalPosterior" << count << ".jpg";
		imwrite( save_totalPosterior.str(), postThreshold);

		//output_cap.write(orig);


		/*
		clock_t tStart = clock();
		cap >> image;
		original = image.clone();
		prob = ProbabilityMap(image, localCues, globalCues, csCues, ApMat);
		frame = BGSubtraction(bg, image, prob);

		imshow( "Display Original", original );
		imshow( "Display BG Mask", frame );
		imshow( "Display Probability", prob );

		totalTime = (double)(clock() - tStart)/CLOCKS_PER_SEC;
		totalTimeInt = totalTime*1000;
		framesSync = totalTimeInt/40;
		wKey = totalTimeInt%40;

		for(int i = 0; i < framesSync; i++)
			cap >> image;

		if(wKey == 0)
			cvWaitKey(1);
		else {
			cvWaitKey(wKey);
		}
		*/
		count = count + 1;
		cvWaitKey(10);

		/*else {
			c = cvWaitKey(wKey);
			if(c==105) {  // ZOOM IN
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&rzoom=1500";
				PTZControl(action);
			}
			if(c==111) {  // ZOOM OUT
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&rzoom=-1500";
				PTZControl(action);
			}
			if(c==114) {  // RIGHT
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&move=right";
				PTZControl(action);
			}
			if(c==108) {  // LEFT
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&move=left";
				PTZControl(action);
			}
			if(c==117) {  // UP
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&move=up";
				PTZControl(action);
			}
			if(c==100) {  // DOWN
				action = "http://root:poseidon@147.27.11.212:80//axis-cgi/com/ptz.cgi?camera=1&move=down";
				PTZControl(action);
			}
		}*/
		bgFlag = 1;
		salientOutput.release();
		frame.release();

	}

	return 0;
}

#endif

