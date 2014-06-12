// Horizon stuff
// 
// Paris Kaimakis Feb 2013



#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "config.h"
#include "foregroundprivates.h"
#include "foreground.h"


using std::cout;
using std::endl;

extern t_settings* g_pSetts;
extern t_horizon* g_pHorizon;


void createHorizon(int U, int V, t_horizon* pHorizon)
{
	pHorizon->affMat = cvCreateMat(2, 3, CV_32FC1) ;
	pHorizon->affMatRot = cvCreateMat(2, 3, CV_32FC1) ;

	pHorizon->ImLandInSea = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);

}


void deleteHorizon(t_horizon* pHorizon)
{
	cvReleaseMat(&pHorizon->affMat);
	cvReleaseMat(&pHorizon->affMatRot);

	cvReleaseImage(&pHorizon->ImLandInSea);
}



void HorizonFindVforThisU(IplImage* Ih, int uCol, int& vColourTransition, int& iTrans)
{

	int V = Ih->height;
	int U = Ih->width;
	int ws = Ih->widthStep;
	
	unsigned char* aRowIh  = (unsigned char*) Ih->imageData + ws;
	int colourRightPrev = 0;
	vColourTransition = V-1;

	int v0colour = (int) aRowIh[uCol];
	iTrans = v0colour;

	for(int v=0; v<V; ++v)
	{
		int colour =  (int) aRowIh[uCol];

		if(colour!=v0colour)
		{
			vColourTransition = v;
			return;
		}
			
		aRowIh += ws;
	}

}



int fitLineWithRegression(int n, double ax[], double ay[], int abValid[], double& alpha, double& beta, double& errorRMS)
{
	// There's n entries in the arrays. Entries whose abValid value is 0
	//	are ommitted.


	double Sx=0, Sy=0;
	double Sxx=0, Sxy=0, Syy=0;
	int nReal = 0;
	for(int i=0; i<n; ++i)
	{
		if(abValid[i])
		{
			++nReal;
			double x = ax[i];
			double y = ay[i];
			
			Sx += x;
			Sy += y;

			Sxx += (x*x);
			Sxy += (x*y);
			Syy += (y*y);
		}
	}


	if(nReal<2)
	{
		//cout << "Insufficient number of valid points found." << endl;
		errorRMS = LARGE_VALUE;
		return 0;
	}



	beta = (nReal*Sxy - Sx*Sy)/(nReal*Sxx-Sx*Sx);
	alpha = Sy/nReal - beta*Sx/nReal;


	// Calculate sqErr:
	double sqErr=0;
	for(int i=0; i<n; ++i)
	{
		if(abValid[i])
		{
			double err = ay[i] - (alpha + beta*ax[i]);
			sqErr += (err*err);

		}
	}

	errorRMS = sqrt(sqErr);



	return 1;
}





double locateHorizonEndpoints(IplImage* Ih, double& u1, double& v1, double& u2, double& v2)
{
	int U = Ih->width;
	int V = Ih->height;

	int n = 10;
	int nPxlOffset = 5;
	double* av = (double*) calloc(n, sizeof(double));
	double* au = (double*) calloc(n, sizeof(double));
	int* abValid = (int*) calloc(n, sizeof(int));
	int* aTopColour = (int*) calloc(n, sizeof(int));	// ==1 if white-to-black, ==-1 black-to-white
	int Du = ((double) U-1-2*nPxlOffset)/(n-1);


	int nTopWhite=0, nTopBlack=0;

	for(int i=0; i<n; ++i)
	{
		abValid[i] = 1;
		au[i] = nPxlOffset + i*Du;
		int iv;
		HorizonFindVforThisU(Ih, au[i], iv, aTopColour[i]);
		av[i] = (double) iv;

		if(iv==0 || iv==V-1)
			abValid[i] = 0;

		

		/*
		cout << "u[" << i << "]=" << au[i] << "\t";
		cout << "\tv[" << i << "]=" << av[i] << "\t";
		cout << "\tbValid[" << i << "]=" << abValid[i] << "\t";
		//cout << "\tiTopCol[" << i << "]=" << aTopColour[i] << "\t";
		cout << endl;
		//*/

		// Keep counts:
		if(aTopColour[i] == 255)
			++nTopWhite;
		if(aTopColour[i] == 0)
			++nTopBlack;

	}



	for(int i=0; i<n; ++i)
	{
		if(nTopWhite>=nTopBlack)
		{
			if(aTopColour[i]==0)
				abValid[i]=0;
		}
		else
		{
			if(aTopColour[i]==255)
				abValid[i]=0;
		}


		//cout << "\tbValidAftr[" << i << "]=" << abValid[i] << "\t";

	}


	// Linear Regression:
	double alpha, beta, errorRMS;
	if(!fitLineWithRegression(n, au, av, abValid, alpha, beta, errorRMS))
	{
		alpha = V-1;
		beta = 0;
	}



	u1 = 0;
	v1 = alpha + beta*u1;

	u2 = U-1;
	v2 = alpha + beta*u2;


	/*
	// Draw?
	// Draw the points:
	for(int i=0; i<n; ++i)
		if(abValid[i])
			cvCircle(Ih, cvPoint(au[i], av[i]), 2, cvScalar(100));

	// Draw Line:
	cvLine(Ih, cvPoint(u1,v1), cvPoint(u2,v2), cvScalar(180), 1);
	cvShowImage("Regr", Ih);
	//*/


	free(av);
	free(au);
	free(abValid);

	return errorRMS;
}



void makeGaussianTransform(IplImage* Isrc, IplImage* Idst)
{
	// same size, 1chan, 32bit

	int U = Isrc->width;
	int V = Isrc->height;
	int ws = Isrc->widthStep;
	float* aRowIsrc = (float*) Isrc->imageData;
	float* aRowIdst = (float*) Idst->imageData;
	float sigma = .03;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			float x = aRowIsrc[u];
			aRowIdst[u] = exp(-x*x/(2*sigma*sigma));

		}

		aRowIsrc += U;
		aRowIdst += U;
	}


}




void trainGaussian(IplImage* It, gsl_vector* m, gsl_matrix* Sinv)
{
	// It is a template image. Estimate the distribution of pixel colours.
	// It is 8bit 3chan

	int U = It->width;
	int V = It->height;
	int nPxls = U*V;
	int ws = It->widthStep;
	unsigned char* aRowI = (unsigned char*) It->imageData;

	// Todo:
	gsl_vector* x = gsl_vector_calloc(3);
	gsl_matrix* S = gsl_matrix_calloc(3,3);
	gsl_matrix* eMat = gsl_matrix_calloc(3,1);
	gsl_matrix* ee = gsl_matrix_calloc(3,3);
	gsl_matrix* SLU = gsl_matrix_calloc(3,3);
	gsl_permutation* p3 = gsl_permutation_calloc(3);


	gsl_vector_set_zero(m);
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			int uTimes3 = u*3;
			double B = (double) aRowI[uTimes3+0];
			double G = (double) aRowI[uTimes3+1];
			double R = (double) aRowI[uTimes3+2];

			x->data[0] = B;
			x->data[1] = G;
			x->data[2] = R;

			gsl_vector_add(m, x);

		}

		aRowI += ws;
	}
	gsl_vector_scale(m, 1.0/nPxls);




	// 1) Covariance:
	gsl_matrix_set_zero(S);
	gsl_matrix_set_zero(Sinv);
	aRowI = (unsigned char*) It->imageData;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			int uTimes3 = u*3;
			double B = (double) aRowI[uTimes3+0];
			double G = (double) aRowI[uTimes3+1];
			double R = (double) aRowI[uTimes3+2];

			
			eMat->data[0] = B - m->data[0];
			eMat->data[1] = G - m->data[1];
			eMat->data[2] = R - m->data[2];

			gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 , eMat, eMat, 0.0, ee);
			gsl_matrix_add(S, ee);

		}

		aRowI += ws;
	}
	gsl_matrix_scale(S, 1.0/(nPxls-1));

	
	S->data[1] = 0;
	S->data[2] = 0;
	S->data[3] = 0;
	S->data[5] = 0;
	S->data[6] = 0;
	S->data[7] = 0;	
	S->data[0] = 1000;
	S->data[4] = 1000;
	S->data[8] = 1000;


	int s;
	gsl_matrix_memcpy(SLU, S);
	gsl_linalg_LU_decomp(SLU, p3, &s);
	double Sdet = gsl_linalg_LU_det(SLU, s);
	gsl_linalg_LU_invert(SLU, p3, Sinv);






	gsl_vector_free(x);
	gsl_matrix_free(S);
	gsl_matrix_free(eMat);
	gsl_matrix_free(ee);
	gsl_matrix_free(SLU);
	gsl_permutation_free(p3);
}




void getLikelihood(IplImage* Ibgr, gsl_vector* m, gsl_matrix* Sinv, IplImage* Ilike)
{

	int U = Ibgr->width;
	int V = Ibgr->height;
	int wsIbgr = Ibgr->widthStep;
	int nPxls = U*V;
	unsigned char* aRowI = (unsigned char*) Ibgr->imageData;
	float* aRowIlike = (float*) Ilike->imageData;


	t_nrml gaussian;
	gaussian.m = m;
	gaussian.Sinv = Sinv;



	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{	
			int uTimes3 = u*3;
			double B = (double) aRowI[uTimes3++];
			double G = (double) aRowI[uTimes3++];
			double R = (double) aRowI[uTimes3];

			aRowIlike[u] = pixelColourLikelihoodFast(B, G, R, &gaussian);

		}

		aRowI += wsIbgr;
		aRowIlike += U;

	}



}




void makeHorizonHorizontalPaddedWay(t_horizon* pHorizon, IplImage* I, IplImage* Ipadded, IplImage* Irot)
{
	int xBoxLeft = pHorizon->xBoxLeft;
	int yBoxTop  = pHorizon->yBoxTop;
	int U = pHorizon->U;
	int V = pHorizon->V;
	CvMat* affMatNew = pHorizon->affMatRot;

	cvSet(Ipadded, cvScalarAll(255));
	CvRect centralRect;
	centralRect = cvRect(xBoxLeft, yBoxTop, U, V);
	cvSetImageROI(Ipadded, centralRect);
	cvCopy(I, Ipadded);
	cvResetImageROI(Ipadded);


	cvWarpAffine(Ipadded, Irot, affMatNew, CV_WARP_FILL_OUTLIERS, cvScalarAll(255));
	
}





void findHorizon(CvCapture* pCapt, IplImage* Iraw, IplImage* I, t_horizon* pHorizon)
{
	// 1) Likelihoods
	// 2) Posteriors
	// 3) Smooths
	// 4) Thresholds, Line-fit
	// 5) Choose the final horizon
	// 6) Angle and distance from top of screen
	// 7) Affine matrix forward + inverse rotation	
	// 8) Calculate Mask for Land-in-sea (ImLandInSea)
	// 9) U, V, xM, yM, affMat for new (larger) rotated image

	
	cout << "\tLocating horizon... " ;
	cout.flush();


	double& xM = pHorizon->xM;
	double& yM = pHorizon->yM;
	double& thetaDeg = pHorizon->thetaHdeg;
	CvMat* affMat = pHorizon->affMat;
	IplImage* ImLandInSea = pHorizon->ImLandInSea;


	// Usefuls:
	int U = I->width;
	int V = I->height;
	int h = (int) (V/10 + .5);
	int hx = h;
	int hy = h;
	CvSize size = cvSize(U,V);

	IplImage* Itemplate		= cvCreateImage(cvSize(hx, hy), IPL_DEPTH_32F, 3);
	IplImage* Itempl8U		= cvCreateImage(cvSize(hx, hy), IPL_DEPTH_8U, 3);
	IplImage* ItmRes2small	= cvCreateImage(cvSize(U-hx+1,V-hy+1), IPL_DEPTH_32F, 1);
	IplImage* I1chan		= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* ItmRes8U		= cvCreateImage(size, IPL_DEPTH_8U, 1); 
	IplImage* IlikeSea		= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IlikeSky		= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* Ipost			= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* Isum			= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostSky		= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostSea		= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostSeaCpy	= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostThresh	= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostSeaRot32F= cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* IpostSeaRot	= cvCreateImage(size, IPL_DEPTH_8U, 1);
	IplImage* ImSkyFromHorizon = cvCreateImage(size, IPL_DEPTH_8U, 1);
	IplImage* ImSkyFromHorizonLarger = cvCreateImage(cvSize(U+2, V+2), IPL_DEPTH_8U, 1);




	gsl_vector* mSky = gsl_vector_calloc(3);
	gsl_matrix* SinvSky = gsl_matrix_calloc(3,3);	
	gsl_vector* mSea = gsl_vector_calloc(3);
	gsl_matrix* SinvSea = gsl_matrix_calloc(3,3);	


	double scFctr = (double)V/180;



	// 0) Capture a new (nonrotated) frame:	//
	Iraw = cvQueryFrame(pCapt);
	int Ubig = Iraw->width;
	int Vbig = Iraw->height;
	cvSetImageROI(Iraw, cvRect(1, 1, Ubig-2, Vbig-2));
	cvResize(Iraw, I);
	cvResetImageROI(Iraw);



	// 1) Likelihoods			// 

	int top = (int) (.0*V);
	int bottom = top + hy;
	int left = (int) (.8*U);
	int right = left + hx;
	CvPoint q1 = cvPoint(left, top);
	CvPoint q2 = cvPoint(right, bottom);
	cvSetImageROI(I, cvRect(left, top, hx, hy));
	cvCopy(I, Itempl8U);
	cvResetImageROI(I);
	trainGaussian(Itempl8U, mSky, SinvSky);
	getLikelihood(I, mSky, SinvSky, IlikeSky);


	CvPoint p1 = cvPoint(.05*U-hx/2, .85*V); // left-centre
	CvPoint p2 = cvPoint(p1.x+hx, p1.y+hy);
	cvSetImageROI(I, cvRect(p1.x, p1.y, hx, hy));
	cvCopy(I, Itempl8U);
	cvResetImageROI(I);
	trainGaussian(Itempl8U, mSea, SinvSea);
	getLikelihood(I, mSea, SinvSea, IlikeSea);



	normalise32FImage(IlikeSea, 0.00, 1);
	normalise32FImage(IlikeSky, 0.01, 1);
	cvAdd(IlikeSea, IlikeSky, Isum);
	cvDiv(IlikeSea, Isum, IpostSea);
	//cvShowImage("2) Posterior Sea",  IpostSea);
	cvCopy(IpostSea, IpostSeaCpy);



	normalise32FImage(IlikeSea, 0.01, 1);
	normalise32FImage(IlikeSky, 0.00, 1);
	cvAdd(IlikeSea, IlikeSky, Isum);
	cvDiv(IlikeSky, Isum, IpostSky);




	cvSetImageROI(IpostSky, cvRect(1, 1, U-2, V-2));
	cvSmooth(IpostSky, IpostSky, CV_GAUSSIAN, 21, 21, 10.0, 10.0);
	cvResetImageROI(IpostSky);

	cvSetImageROI(IpostSea, cvRect(1, 1, U-2, V-2));
	cvSmooth(IpostSea, IpostSea, CV_GAUSSIAN, 21, 21, 10.0, 10.0);
	cvResetImageROI(IpostSea);



	// 4) Threshold, Line Fit :	
	double x1=0, x2=U-1;

	double x1sky=0, y1sky=0, x2sky=U-1, y2sky=0;
	cvThreshold(IpostSky, IpostThresh, 0.5, 1.0, CV_THRESH_BINARY); // post 0.5, likelihood 0.9
	cvConvertScale(IpostThresh , ItmRes8U, 255.0);	
	//cvShowImage("Thresh Sky", ItmRes8U);
	double errSkyRMS = locateHorizonEndpoints(ItmRes8U, x1, y1sky, x2, y2sky);
	cvLine(I, cvPoint(x1sky, y1sky), cvPoint(x2sky, y2sky), cvScalar(255, 100, 100), 2);


	double x1sea=0, y1sea=0, x2sea=U-1, y2sea=0;
	cvThreshold(IpostSea, IpostThresh, 0.5, 1.0, CV_THRESH_BINARY); // post 0.5, likelihood 0.9
	cvConvertScale(IpostThresh , ItmRes8U, 255.0);	
	//cvShowImage("Thresh Sea", ItmRes8U);
	double errSeaRMS = locateHorizonEndpoints(ItmRes8U, x1, y1sea, x2, y2sea);
	cvLine(I, cvPoint(x1sea, y1sea), cvPoint(x2sea, y2sea), cvScalar(255, 0, 0), 2);

	//cout << "errSkyRMS=" << errSkyRMS << "	errSeaRMS=" << errSeaRMS << endl;


	//------------------------------//
	double y1, y2;
	double errRMS=0;
	int thErr = 15*scFctr;
	int bHorizonLocated = 0;
	if(errSkyRMS<thErr || errSeaRMS<thErr)
	{
		// Horizon fit is acceptable, choose the best:
		if(errSkyRMS<errSeaRMS)
		{
			errRMS = errSkyRMS;
			y1 = y1sky;
			y2 = y2sky;
		}
		else
		{
			errRMS = errSeaRMS;
			y1 = y1sea;
			y2 = y2sea;
		}

		bHorizonLocated = 1;
	}
	else
	{
		cout << "\t\t\tBad horizon fit. Switching to non-horizon mode..." << endl;
		y1 = V-1;
		y2 = V-1;
		bHorizonLocated = 0;
		g_pSetts->sigmaMin = 3.0;
	}
	
	cvLine(I, cvPoint(x1, y1), cvPoint(x2, y2), cvScalar(0, 100, 255), 2);



	double Dy = y2-y1;
	double Dx = x2-x1;
	
	xM = (x1 + x2)/2;
	yM = (y1 + y2)/2;

	


	thetaDeg = 180/PI*atan(Dy/Dx);

	//cvCircle(I, cvPoint(xM,yM), 3, cvScalar(0,255,0), -1);
	cvRectangle(I, p1, p2, cvScalar(0, 255, 255), 1);
	cvRectangle(I, q1, q2, cvScalar(0, 255, 0), 1);
	//cvShowImage("2) Horizon Analysis", I);
	//*/



	// 7) Affine matrix forward + inverse rotation	//

	CvPoint2D32f p;
	p.x = xM;
	p.y = yM;


	cv2DRotationMatrix(p, thetaDeg, 1.0, affMat); 



	cvSet(ImLandInSea, cvScalarAll(0));
	if(bHorizonLocated)
	{
		cvSet(ImSkyFromHorizon, cvScalarAll(0));
		cvSet(ImSkyFromHorizonLarger, cvScalarAll(0));
		int pushHorDown = (int) std::max(errRMS +0.5, 5*scFctr);

		cvLine(ImSkyFromHorizon, cvPoint(x1, y1+pushHorDown ), cvPoint(x2, y2+pushHorDown ), cvScalar(255), 2);
		int flagsForFloodFill = 4| CV_FLOODFILL_MASK_ONLY| CV_FLOODFILL_FIXED_RANGE | (255<<8);
		int yMforFloodfill = (int) std::max(0.0, (double) (yM-3));
		cvFloodFill(ImSkyFromHorizon, cvPoint(xM, yMforFloodfill), cvScalarAll(255), cvScalarAll(0), cvScalarAll(0), NULL, flagsForFloodFill, ImSkyFromHorizonLarger);
		cvSetImageROI(ImSkyFromHorizonLarger, cvRect(1,1,U, V));

		cvSet(IpostSeaCpy, cvScalar(1), ImSkyFromHorizonLarger);
		cvLine(IpostSeaCpy, cvPoint(x1, y1+pushHorDown ), cvPoint(x2, y2+pushHorDown ), cvScalar(1), 3);

	

		cvCopy(IpostSeaCpy, IpostSeaRot32F);
		cvThreshold(IpostSeaRot32F, IpostSeaCpy, 0.5, 1.0, CV_THRESH_BINARY);
		cvScale(IpostSeaCpy, IpostSeaCpy, -1.0, 1.0);
		cvConvertScale(IpostSeaCpy, IpostSeaRot, 255.0);	
		cvDilate(IpostSeaRot, IpostSeaRot);
		//cvShowImage("2 Dilated", IpostSeaRot);
		cvErode(IpostSeaRot, IpostSeaRot);
		//cvShowImage("3 Eroded", IpostSeaRot);
		cvDilate(IpostSeaRot, IpostSeaRot);
		//cvDilate(IpostSeaRot, IpostSeaRot);
		//cvShowImage("Sea rot", IpostSeaRot);
		
		cvCopy(IpostSeaRot, ImLandInSea);

	}


	//----------------------------------------------------------//

	double cosTheta = cos(fabs(thetaDeg)*PI/180);
	double sinTheta = sin(fabs(thetaDeg)*PI/180);
	int Vrot = (int) (V*cosTheta + U*sinTheta + .5);
	int Urot = (int) (U*cosTheta + V*sinTheta + .5);


	double Dx1, Dy1, alpha, xi, d, xMdash;
	double Dx2, Dy2, alpha2, zeta, d2, yMdash;


	if(thetaDeg>0)
	{
		Dx1 = xM-0;
		Dy1 = yM-0;

		// yMdash governed by top-right
		Dx2 = U-xM;
		Dy2 = yM-0;
	}
	else
	{
		Dx1 = xM-0;
		Dy1 = V-yM;

		Dx2 = xM-0;
		Dy2 = yM-0;
	}

	double thetaDegFabs = fabs(thetaDeg);

	alpha = atan(Dy1/Dx1)*180/PI;
	xi = fabs(thetaDegFabs - alpha);
	d = sqrt(Dx1*Dx1 + Dy1*Dy1);
	xMdash = d*cos(xi*PI/180);

	alpha2 = atan(Dy2/Dx2)*180/PI;
	zeta = fabs(90 - alpha2 - thetaDegFabs);
	d2 = sqrt(Dx2*Dx2 + Dy2*Dy2);
	yMdash = d2*cos(zeta*PI/180);
	

	int xBoxLeft = (int) (xMdash - xM + .5);
	int yBoxTop =  (int) (yMdash - yM + .5);

	if(xBoxLeft<0)
	{
		Urot -= (--xBoxLeft); // (increment) the -- is related to the +.5 which was added as rounding assuming positive value
		xBoxLeft = 0;
		xMdash = xM;
	}
	if(yBoxTop<0)
	{
		Vrot -= (--yBoxTop); // (increment) 
		yBoxTop = 0;
		yMdash = yM;
	}

	CvMat* affMatNew = cvCreateMat(2, 3, CV_32FC1);
	CvPoint2D32f pNew;
	pNew.x = xMdash;
	pNew.y = yMdash;
	cv2DRotationMatrix(pNew, thetaDeg, 1.0, affMatNew); 

	pHorizon->U = U;
	pHorizon->V = V;
	pHorizon->Urot = Urot;
	pHorizon->Vrot = Vrot;
	pHorizon->xMrot = xMdash;
	pHorizon->yMrot = yMdash;
	pHorizon->xBoxLeft = xBoxLeft;
	pHorizon->yBoxTop = yBoxTop;
	cvCopy(affMatNew, pHorizon->affMatRot);

	
	IplImage* Ipadded = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 3);
	IplImage* Irot = cvCreateImage(cvSize(Urot, Vrot), IPL_DEPTH_8U, 3);
	cvSetZero(Ipadded);
	makeHorizonHorizontalPaddedWay(pHorizon, I, Ipadded, Irot);



	cvReleaseImage(&Irot);
	cvReleaseImage(&Ipadded);
	cvReleaseMat(&affMatNew);
	


	cvWaitKey(1);




	// Releases:
	cvReleaseImage(&IlikeSea);
	cvReleaseImage(&IlikeSky);
	cvReleaseImage(&Ipost);
	cvReleaseImage(&Isum);
	cvReleaseImage(&IpostThresh);
	cvReleaseImage(&IpostSky);
	cvReleaseImage(&IpostSea);
	cvReleaseImage(&IpostSeaCpy);
	cvReleaseImage(&IpostSeaRot32F);
	cvReleaseImage(&IpostSeaRot);
	cvReleaseImage(&ImSkyFromHorizon);
	cvReleaseImage(&ImSkyFromHorizonLarger);

	gsl_vector_free(mSky);
	gsl_vector_free(mSea);
	gsl_matrix_free(SinvSky);
	gsl_matrix_free(SinvSea);


	if(bHorizonLocated)
		cout << "Success!" << endl;
	else
		cout << "Failed!" << endl;


}





void makeHorizonHorizontal(IplImage* Isrc, CvMat* affMat, IplImage* Idst)
{

	int U = Isrc->width;
	int V = Isrc->height;
	

	cvWarpAffine(Isrc, Idst, affMat,CV_WARP_FILL_OUTLIERS, cvScalarAll(1));

	//cvLine(I1, cvPoint(0, d), cvPoint(U-1, d), cvScalar(255, 255,0), 1);
	//cvShowImage("Rotated", Idst);


}





void displayUnrotatedImage(char* szText, IplImage* Irot)
{
	t_horizon* pHorizon = g_pHorizon;

	IplImage* Iunrot = cvCloneImage(Irot);//CreateImage(cvSize(320,180), Irot->depth, Irot->nChannels);

	cvWarpAffine(Irot, Iunrot, pHorizon->affMatRot, CV_WARP_INVERSE_MAP|CV_WARP_FILL_OUTLIERS, cvScalarAll(0));

	cvSetImageROI(Iunrot, cvRect(pHorizon->xBoxLeft, pHorizon->yBoxTop, pHorizon->U, pHorizon->V));

	cvShowImage(szText, Iunrot);

	cvReleaseImage(&Iunrot);
}



void returnUnrotatedImage(char* szText, IplImage* Irot, IplImage* Iunrot)
{
	t_horizon* pHorizon = g_pHorizon;

	IplImage* IunrotTemp = cvCloneImage(Irot);//CreateImage(cvSize(320,180), Irot->depth, Irot->nChannels);

	cvWarpAffine(Irot, IunrotTemp, pHorizon->affMatRot, CV_WARP_INVERSE_MAP|CV_WARP_FILL_OUTLIERS, cvScalarAll(0));

	cvSetImageROI(IunrotTemp, cvRect(pHorizon->xBoxLeft, pHorizon->yBoxTop, pHorizon->U, pHorizon->V));

	cvCopy(IunrotTemp, Iunrot);

	cvShowImage(szText, IunrotTemp);
	cvShowImage("Iunrot", Iunrot);


	cvReleaseImage(&IunrotTemp);
}


