// Manage and display the track of the targets in the scene.
// Paris Kaimakis 23 May 2013


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "foreground.h"
#include "foregroundprivates.h"
#include "foregroundconstants.h"



using std::cout;
using std::endl;



#ifdef F_TRACK

void createTrackArray(CvCapture* pCapture, int U, int V, t_track aTrack[])
{
	int duration = cvGetCaptureProperty(pCapture, CV_CAP_PROP_FRAME_COUNT);


	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		t_track* pTrack = &aTrack[iT];

		pTrack->aU = (double*) std::calloc(duration, sizeof(double));
		pTrack->aV = (double*) std::calloc(duration, sizeof(double));

		pTrack->Itrack = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 3);
		cvSet(pTrack->Itrack, cvScalarAll(255));
	}


}



void deleteTrackArray(t_track aTrack[])
{
	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		t_track* pTrack = &aTrack[iT];

		free(pTrack->aU);
		free(pTrack->aV);
		cvReleaseImage(&pTrack->Itrack);
	}
}




void drawTrackArrayFrameByFrame(int iFrame, t_track aTrack[])
{
	// I is 3chan 8bit.


	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		t_track* pTrack = &aTrack[iT];
		IplImage* Itrack = pTrack->Itrack;

		std::stringstream ssNum;
		ssNum << iT;
		std::string szWinName = "Track " + ssNum.str() ;


		// Too soon to start drawing?
		if(iFrame<1) return;

		double uNow = pTrack->aU[iFrame];
		double vNow = pTrack->aV[iFrame];
		double uBef = pTrack->aU[iFrame-1];
		double vBef = pTrack->aV[iFrame-1];

		// No target now or before?
		if( (uNow==0 && vNow==0) || (uBef==0 && vBef==0) ) return;

		double Du = uNow - uBef;
		double Dv = vNow - vBef;
		double lengthSq = Du*Du+Dv*Dv;
		if(lengthSq>9) return;

		cvLine(Itrack, cvPoint(uBef, vBef), cvPoint(uNow, vNow), cvScalar(0, 0, 0), 3);


		cvShowImage(szWinName.c_str(), Itrack);

	}
}



void drawTrackArrayOneOff(IplImage* Ibgr, int nFrames, t_track aTrack[])
{


	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		t_track* pTrack = &aTrack[iT];
		IplImage* Itrack = cvCloneImage(Ibgr);

		std::stringstream ssNum;
		ssNum << iT;
		std::string szWinName = "One-Off Track " + ssNum.str() ;

		for(int iFrame=268; iFrame<nFrames; ++iFrame)
		{
			if(iFrame>=1)
			{

				double uNow = pTrack->aU[iFrame];
				double vNow = pTrack->aV[iFrame];
				double uBef = pTrack->aU[iFrame-1];
				double vBef = pTrack->aV[iFrame-1];



				// Target is present
				if( (uNow!=0 && vNow!=0) && (uBef!=0 && vBef!=0) ) 
				{

					double Du = uNow - uBef;
					double Dv = vNow - vBef;
					double lengthSq = Du*Du+Dv*Dv;
					if(lengthSq<9)
					{

						cvLine(Itrack, cvPoint(uBef, vBef), cvPoint(uNow, vNow), cvScalar(0, 0, 0), 2);
					
						cvShowImage(szWinName.c_str(), Itrack);
						cvWaitKey(1);
					}
				}
			}
		}


		cvShowImage(szWinName.c_str(), Itrack);

		cvReleaseImage(&Itrack);

	}
}



void tidyTrackArray(int nTargets, int iFrame, t_track aTrack[])
{

	if(iFrame<1) return;

	double distThresh = ((double )aTrack[0].Itrack->height)/20;
	double distSqThresh = 9;//distThresh*distThresh;


	int aTargetAssociation[MAX_TARGETS] = {0};
	for(int iT=0; iT<MAX_TARGETS; ++iT)
		aTargetAssociation[iT] = -1;


	for(int iTnow=0; iTnow<nTargets; ++iTnow)
	{

		double uNow = aTrack[iTnow].aU[iFrame];
		double vNow = aTrack[iTnow].aV[iFrame];


		double minDistSq = LARGE_VALUE;
		for(int iTbef=0; iTbef<MAX_TARGETS; ++iTbef)
		{

			if(aTrack[iTbef].aU[iFrame-1] == 0) break;

			double uBef = aTrack[iTbef].aU[iFrame-1];
			double vBef = aTrack[iTbef].aV[iFrame-1];

			double Du = uNow-uBef;
			double Dv = vNow-vBef;

			double distSq = Du*Du + Dv*Dv;

			if(distSq<minDistSq && distSq<distSqThresh)
			{
				minDistSq = distSq;
				//cout << sqrt(minDistSq) << endl;
				aTargetAssociation[iTnow] = iTbef;
			}


		}

	}



	double aUtmp[MAX_TARGETS];
	double aVtmp[MAX_TARGETS];
	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		aUtmp[iT] = aTrack[iT].aU[iFrame];
		aVtmp[iT] = aTrack[iT].aV[iFrame];
	}


	// Switch:
	for(int iT=0; iT<MAX_TARGETS; ++iT)
	{
		int iTtoSwitchTo = aTargetAssociation[iT];
	
		if(iTtoSwitchTo >= 0)
		{
			aTrack[iT].aU[iFrame] = aUtmp[iTtoSwitchTo];
			aTrack[iT].aV[iFrame] = aVtmp[iTtoSwitchTo];
		}
	}

}

#endif
