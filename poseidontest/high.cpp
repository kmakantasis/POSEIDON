// Poseidon main function
//
// Paris Kaimakis 12 Oct 2012

#include <iostream>
#include <fstream>
#include <ctime>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>
#include "foreground.h"
#include "datasetinfo.h"
#include "foregroundinterface.h"
#include "videocapture.h"
#include "poseidontest.h"


using std::cout;
using std::endl;

// Switch one of the following on:
//#define LF_FULL_COVAR_TOO
#define LF_DIAG_COVAR_TOO
//#define F_PAPER

extern t_settings* g_pSetts;



int doEverything(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[])
{

	//if(!getNewFrame(pVid))
		//return 0;




	int& iFrame = pVid->iFrame;
	++iFrame;
	//cout << "Frame " << iFrame << endl;



	//if(iFrame<1100 || iFrame % 1 !=0) return 1;


	//cout << "nframes = " << cvGetCaptureProperty(pVid->capture, CV_CAP_PROP_FRAME_COUNT) << endl;

#ifdef F_PAPER
	displayPixelIntensityGraphsForPaper(pVid->Iresults, pVid->iFrame, &pBg->paper);
#endif




	makeHorizonHorizontalPaddedWay(pHorizon, pVid->Iresized, pBg->Iresized, pBg->ItoProcess);
	cvCopy(pBg->Iresized, pBg->Iresults);
	cvCopy(pBg->ItoProcess, pBg->Iresized);




	showBackgroundAsImage(pBg);

	
#ifdef LF_FULL_COVAR_TOO
	removeBackground(&pBg->mons, pBg->ImLandInSea, pHorizon->yM, pBg->aRGB, pBg->ItoProcess, pBg->Ilhood, pBg->ImaskSky, pBg->IforSegmo, pBg->IjointBG);
#endif
#ifdef LF_DIAG_COVAR_TOO
	removeBackground(&pBg->mons, pBg->ImLandInSea, pHorizon->yMrot, pBg->aRGBdiag, pBg->ItoProcess, pBg->IlhoodDiag, pBg->ImaskSky, pBg->IforSegmo, pBg->IjointBG);
#endif


		
	
#ifdef LF_FULL_COVAR_TOO
	int nTargets = postprocess(pHorizon, pBg->Ilhood, pBg->ImLandInSea, pVid->ItoProcess, &pBg->mons, pBg->kernel, pBg->IforSegmo, pVid->Iresults, pBg->IsSignificant, iFrame, aTrack);
#endif
#ifdef LF_DIAG_COVAR_TOO
	int nTargets = postprocess(pHorizon, pBg->IlhoodDiag, pBg->ImLandInSea, pBg->ItoProcess, &pBg->mons, pBg->kernel, pBg->IforSegmo, pBg->Iresults, pBg->IsSignificant, iFrame, aTrack);
#endif
	//*/

	//if(iFrame%5 == 0)
	adaptBackground(&pBg->mons, 0.01, pBg->ItoProcess, pBg);

	/*
	std::ofstream fout;
	fout.open("Targets.txt",std::ios_base::app);
	fout << "" << endl; 
	
	fout << "Frame=" << iFrame << "	nT=" << nTargets;
	fout.close();
	//*/

	//if(nTargets>0)
		//cvWriteFrame(pVid->writer, pVid->Iresized);
	//cvWriteFrame(pVid->writer, pVid->Iresults);

	//if(iFrame==401||iFrame==875||iFrame==1137||iFrame==1326||iFrame==1757||iFrame==1966||iFrame==2520||iFrame==3150||iFrame==4247||iFrame==4473||iFrame==5098||iFrame==5358)
	//if(iFrame==249)
		//cvWaitKey();

	
	// Draw target track:
	//if(iFrame>=268)
		//drawTrackArrayFrameByFrame(iFrame, aTrack);
	
	/*
	if(iFrame==459)	
	{
		drawTrackArrayOneOff(pBg->Iresized, pVid->iFrame, aTrack);
		cvWaitKey();
	}
	//*/

	int bContinue = handleKeyInput(pInfo, pVid, pHorizon, pBg);


	return bContinue;
}





void processVideo(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[])
{


	loadVideo(pVid);
	findHorizon(pVid->capture, pVid->Iframe, pVid->Iresized, pHorizon);

	//createBg(pBg, g_pSetts->U, g_pSetts->V);
	createBg(pBg, pHorizon->Urot, pHorizon->Vrot);



	// Temp:
	int Urot = pHorizon->Urot;
	int Vrot = pHorizon->Vrot;
	IplImage* ImPadded = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 1);
	IplImage* ImLandInSeaRot = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 1);
	//cvSet(ImPadded, cvScalarAll(255));
	makeHorizonHorizontalPaddedWay(pHorizon, pHorizon->ImLandInSea, ImPadded, pBg->ImLandInSea);
	//cvShowImage("ImLandInSeaRot", pBg->ImLandInSea);
	//cvShowImage("ImLandInSea", pHorizon->ImLandInSea);
	//cvWaitKey();


	trainBgModel(pVid->capture, pVid->Iframe, pHorizon, pVid->Iresized, pHorizon->affMat, pBg->ItoProcess, pBg);



	clock_t t1 = clock();
	// Main Loop
	cout << "\tLocating targets... " ;
	cout.flush();

	int bContinue = 1;
	while(bContinue) {
		bContinue = doEverything(pInfo, pVid, pHorizon, pBg, aTrack);

	}
	cout << "Done!" << endl;
	clock_t t2 = clock();
	//cout << "Time elapsed = " << ((float)(t2 - t1))/CLOCKS_PER_SEC << endl;

		
	//drawTrackArrayOneOff(pBg->Iresized, pVid->iFrame, aTrack);

	//destroyBg(pBg);
	//deleteTrackArray(aTrack);


	//unloadVideo(pVid);




}




void initialiseEverythingParis(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[])
{
	// Initialises everything for the sake of integration

	loadVideo(pVid);
	findHorizon(pVid->capture, pVid->Iframe, pVid->Iresized, pHorizon);

	//createBg(pBg, g_pSetts->U, g_pSetts->V);
	createBg(pBg, pHorizon->Urot, pHorizon->Vrot);



	// Temp:
	int Urot = pHorizon->Urot;
	int Vrot = pHorizon->Vrot;
	IplImage* ImPadded = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 1);
	IplImage* ImLandInSeaRot = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 1);
	//cvSet(ImPadded, cvScalarAll(255));
	makeHorizonHorizontalPaddedWay(pHorizon, pHorizon->ImLandInSea, ImPadded, pBg->ImLandInSea);
	//cvShowImage("ImLandInSeaRot", pBg->ImLandInSea);
	//cvShowImage("ImLandInSea", pHorizon->ImLandInSea);
	//cvWaitKey();


	trainBgModel(pVid->capture, pVid->Iframe, pHorizon, pVid->Iresized, pHorizon->affMat, pBg->ItoProcess, pBg);

}






void processVideoTestSet(std::string szPath, int U, int V, t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[])
{

	IplImage* Iuseless = cvCreateImage(cvSize(U,V), IPL_DEPTH_8U, 1);
	char c='0';
	cvShowImage("Targets", Iuseless);	
	
	while(1)
	{
		displayInfo(pInfo);


		int iOption = NUM_DATASETS;
		while(iOption<0 || iOption>=NUM_DATASETS)
		{
			cout << "Select option 0-" << NUM_DATASETS-1 << " or Esc to quit: ";
			c = cvWaitKey();//
			if(c==27) {	cout << endl; break; }
			int iOptionLocal = c - '0';//
			
			if(iOptionLocal>=0 && iOptionLocal<NUM_DATASETS)
			{
				iOption = iOptionLocal;
				cout << iOption << endl;
			}
			else
			{
				cout << c << endl;
			}
		}
		if(c==27) break;


		std::string szFilename = pInfo->aDataset[iOption].szFilename;
		std::string szFilenameWithPath = szPath + szFilename;
		cout << "\nProcessing " << szFilename << " " ;
		cout.flush();
		cout << "(" << pInfo->aDataset[iOption].szDescription << ")" << endl;
		createVideoCapture(szPath, szFilenameWithPath, U, V, pVid);
		//createTrackArray(pVid->capture, U, V, aTrack);



		processVideo(pInfo, pVid, pHorizon, pBg, aTrack);


		deleteVideoCapture(pVid);
		//deleteTrackArray(aTrack);


		cout << "\n\n\n" << endl; 

		iOption=NUM_DATASETS+1;

	} //while 1


	cvReleaseImage(&Iuseless);

}
