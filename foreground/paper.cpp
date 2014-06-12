// Display functions to produce figures for a paper/report.
// Paris Kaimakis 11 May 2013


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "foreground.h"

using std::cout;
using std::endl;


void createStuffForPaper(t_paper* pPaper)
{
	
	pPaper->t = 0;
	pPaper->m = 2; //5
	double m = pPaper->m;



	pPaper->IpixelValuesOverTime = cvCreateImage(cvSize(500, 255*m), IPL_DEPTH_8U, 3);
	cvSet(pPaper->IpixelValuesOverTime , cvScalarAll(255));

	pPaper->IsamplesBG = cvCreateImage(cvSize(255*m, 255*m), IPL_DEPTH_8U, 3);
	pPaper->IsamplesGR = cvCreateImage(cvSize(255*m, 255*m), IPL_DEPTH_8U, 3);
	pPaper->IsamplesRB = cvCreateImage(cvSize(255*m, 255*m), IPL_DEPTH_8U, 3);

	cvSet(pPaper->IsamplesBG, cvScalarAll(255));
	cvSet(pPaper->IsamplesGR, cvScalarAll(255));
	cvSet(pPaper->IsamplesRB, cvScalarAll(255));


}



void deleteStuffForPaper(t_paper* pPaper)
{
	cvReleaseImage(&pPaper->IpixelValuesOverTime);
	cvReleaseImage(&pPaper->IsamplesBG);
	cvReleaseImage(&pPaper->IsamplesGR);
	cvReleaseImage(&pPaper->IsamplesRB);

}



void displayPixelIntensityGraphsForPaper(IplImage* Ibgr, int iFrame, t_paper* pPaper)
{



	IplImage* Idispl = pPaper->IpixelValuesOverTime;
	IplImage* IBG = pPaper->IsamplesBG;
	IplImage* IGR = pPaper->IsamplesGR;
	IplImage* IRB = pPaper->IsamplesRB;
	
	int u = 5; //160
	int v= Ibgr->height-5; //60; //Ibgr->height-5; // 160, 80
	double m = pPaper->m;
	int maxIntensity = 255*m;

	int Dt = 10;
	int& t = pPaper->t;
	t+= Dt;


	// Show the colours so far:
	int B = m*(int) CV_IMAGE_ELEM(Ibgr, unsigned char, v, (u * 3) + 0);
	int G = m*(int) CV_IMAGE_ELEM(Ibgr, unsigned char, v, (u * 3) + 1);
	int R = m*(int) CV_IMAGE_ELEM(Ibgr, unsigned char, v, (u * 3) + 2);
	

	int lineThickness = 4; //4

	//cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Bbefore), cvPoint(t,maxIntensity-B), cvScalar(255, 0, 0), lineThickness);
	cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Bbefore), cvPoint(t,maxIntensity-B), cvScalar(0, 0, 0), lineThickness);
	//cvCircle(Idispl, cvPoint(t, B), 1, cvScalar(255, 0, 0), -1);

	//cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Gbefore), cvPoint(t,maxIntensity-G), cvScalar(0, 255, 0), lineThickness);
	cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Gbefore), cvPoint(t,maxIntensity-G), cvScalar(75, 75, 75), lineThickness);
	//cvCircle(Idispl, cvPoint(t, G), 1, cvScalar(0, 255, 0), -1);
	
	//cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Rbefore), cvPoint(t,maxIntensity-R), cvScalar(0, 0, 255), lineThickness);
	cvLine(Idispl, cvPoint(t-Dt, maxIntensity-pPaper->Rbefore), cvPoint(t,maxIntensity-R), cvScalar(150, 150, 150), lineThickness);
	//cvCircle(Idispl, cvPoint(t, R), 1, cvScalar(0, 0, 255), -1);
	
	
	cvShowImage("p vs t", Idispl);



	int circleThickness = 4; //3


	cvCircle(IBG, cvPoint(B, maxIntensity-G), circleThickness, cvScalar(0,0,255), -1);
	cvCircle(IGR, cvPoint(G, maxIntensity-R), circleThickness, cvScalar(0,0,255), -1);
	cvCircle(IRB, cvPoint(R, maxIntensity-B), circleThickness, cvScalar(0,0,255), -1);
	cvShowImage("B vs G", IBG);
	cvShowImage("G vs R", IGR);
	cvShowImage("R vs B", IRB);
	cvWaitKey(10);


	cvCircle(IBG, cvPoint(B, maxIntensity-G), circleThickness, cvScalar(0), -1);
	cvCircle(IGR, cvPoint(G, maxIntensity-R), circleThickness, cvScalar(0), -1);
	cvCircle(IRB, cvPoint(R, maxIntensity-B), circleThickness, cvScalar(0), -1);
	cvShowImage("B vs G", IBG);
	cvShowImage("G vs R", IGR);
	cvShowImage("R vs B", IRB);





	pPaper->Bbefore = B;
	pPaper->Gbefore = G;
	pPaper->Rbefore = R;




	cvRectangle(Ibgr, cvPoint(u, v), cvPoint(u+1, v+1), cvScalar(100,0,0), 2);
	cvShowImage("A typical bg pixel", Ibgr);





}



