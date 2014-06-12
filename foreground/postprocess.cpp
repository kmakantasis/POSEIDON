// Likelihood image to detection
//
// Paris Kaimakis 4/1/13




#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "config.h"
#include "foreground.h"
#include "foregroundprivates.h"


using std::cout;
using std::endl;


extern t_settings* g_pSetts;


//#define F_DISPLAY_INTERMEDIATES


void keepSignificantClusters(t_monster* pMons, IplImage* Iorig, IplImage* Iseeds, IplImage* Isignificant)
{
	// All are silhouette images of same size.
	// Iseed, Iorig both altered and become useless at the end of this function.


	int U = Iorig->width;
	int V = Iorig->height;
	double minClusterSizeProp = g_pSetts->minClusterSize;
	int Uunrot = g_pSetts->U;
	int Vunrot = g_pSetts->V;
	int nPxlThresh = (int) (Uunrot*Vunrot*minClusterSizeProp *minClusterSizeProp ); 
	int flagsForFloodFill = 4| CV_FLOODFILL_MASK_ONLY| CV_FLOODFILL_FIXED_RANGE | (255<<8);
	CvScalar colBlack = cvScalarAll(0);
	CvScalar colWhite = cvScalarAll(255);
	IplImage* IorigClusterMask = pMons->IorigClusterMask;
	

	cvSet(Isignificant, colWhite);
	unsigned char* aRowSeeds = (unsigned char*) Iseeds->imageData;
	int ws = Iseeds->widthStep;
	CvRect rect = cvRect(1, 1, U, V);


	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRowSeeds[u] == (unsigned char) 0)
			{
				// New cluster spotted.


				cvZero(IorigClusterMask);
				cvFloodFill(Iorig, cvPoint(u, v), colWhite, colBlack, colBlack, NULL, flagsForFloodFill, IorigClusterMask);
				cvSetImageROI(IorigClusterMask, rect);

				cvSet(Isignificant, colBlack, IorigClusterMask);
				cvSet(Iseeds, colWhite, IorigClusterMask);
				cvResetImageROI(IorigClusterMask);


			} // cluster spotted in IseedsCpy
		
		} //u

		aRowSeeds += ws;
	} //v


}



void findBoundingBoxDimensions(IplImage* Imask, int lineThickness, CvPoint& topLeft, CvPoint& bottomRight)
{
	// There's only ONE cluster in Imask.

	int U = Imask->width;
	int V = Imask->height;
	unsigned char* aRowI = (unsigned char*) Imask->imageData;
	int ws = Imask->widthStep;


	int xMin=U, yMin=V;
	int xMax=0, yMax=0;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRowI[u] == (unsigned char) 255)
			{
				xMin = std::min(xMin, u);
				yMin = std::min(yMin, v);
				xMax = std::max(xMax, u);
				yMax = std::max(yMax, v);
			}
		}

		aRowI += ws;
	}

	/*
	xMin = std::max(xMin, 1+lineThickness);
	yMin = std::max(yMin, 1+lineThickness);
	xMax = std::min(xMax, U-2-lineThickness);
	yMax = std::min(yMax, V-2-lineThickness);
	*/
	/*
	cout << "xMin=" << xMin << " xMax=" << xMax << " yMin=" << yMin << " yMax=" << yMax << endl;
	cvShowImage("Mask", Imask);
	cvWaitKey();
	*/

	topLeft = cvPoint(xMin, yMin);
	bottomRight = cvPoint(xMax, yMax);
}


int targetAnalysis(IplImage* Ilhood, IplImage* ImAllegedTarget, IplImage* ImLandInSea, IplImage* Ibgr, CvPoint tl, CvPoint br)
{
	int width = br.x-tl.x;
	int height = br.y-tl.y;
	CvSize size = cvSize(width, height);
	IplImage* IpotentialTarget = cvCreateImage(size, IPL_DEPTH_8U, 3);
	IplImage* IpotentialTargetMask = cvCreateImage(size, IPL_DEPTH_8U, 1);
	IplImage* IsAllegedTarg = cvCreateImage(size, IPL_DEPTH_8U, 1);
	IplImage* IlPatch = cvCreateImage(size, IPL_DEPTH_32F, 1);


	//cvShowImage("Isegm",Ilhood);



	CvRect rect = cvRect(tl.x, tl.y, width, height);
	cvSetImageROI(Ibgr, rect);
	cvCopy(Ibgr, IpotentialTarget);
	cvResetImageROI(Ibgr);
	//cvShowImage("Ipt", IpotentialTarget);


	cvSetImageROI(Ilhood, rect);
	//cvShowImage("Ilhd", Ilhood);
	cvCopy(Ilhood, IlPatch);
	cvResetImageROI(Ilhood);


	cvSetImageROI(ImLandInSea, rect);
	cvCopy(ImLandInSea, IpotentialTargetMask);
	cvResetImageROI(ImLandInSea);
	//cvShowImage("Imask", IpotentialTargetMask);


	cvSetImageROI(ImAllegedTarget, rect);
	cvCopy(ImAllegedTarget, IsAllegedTarg);
	cvResetImageROI(ImAllegedTarget);
	//cvShowImage("Before", IsAllegedTarg);
	maskToSilhouette(IsAllegedTarg);
	//cvShowImage("Isat", IsAllegedTarg);


	
	//cvShowImage("Imask bfr", IpotentialTargetMask);
	cvSet(IpotentialTargetMask, cvScalarAll(255), IsAllegedTarg);
	//cvShowImage("Imask aftr", IpotentialTargetMask);
	//cvWaitKey();


	int U = IpotentialTarget->width;
	int V = IpotentialTarget->height;
	unsigned char* aRowI = (unsigned char*) IpotentialTarget->imageData;
	unsigned char* aRowImask = (unsigned char*) IpotentialTargetMask->imageData;
	float* aRowIlPatch = (float*) IlPatch->imageData;
	int ws = IpotentialTarget->widthStep;

	gsl_vector* m = gsl_vector_calloc(3);
	gsl_vector* x = gsl_vector_calloc(3);
	gsl_matrix* S = gsl_matrix_calloc(3,3);
	gsl_matrix* eMat = gsl_matrix_calloc(3,1);
	gsl_matrix* ee = gsl_matrix_calloc(3,3);
	gsl_matrix* SLU = gsl_matrix_calloc(3,3);
	gsl_permutation* p3 = gsl_permutation_calloc(3);


	// Mean:
	int nPxls = 0;
	double totalWeight = 0;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			//if((int)aRowImask[u] == 0)
			//if(aRowIlPatch[u]>0.5)
			{
				++nPxls;

				int uTimes3 = u*3;
				double B = (double) aRowI[uTimes3+0];
				double G = (double) aRowI[uTimes3+1];
				double R = (double) aRowI[uTimes3+2];

				x->data[0] = B;
				x->data[1] = G;
				x->data[2] = R;

				// scale wrt likelihood:

				float lhoodFg = 1.0-aRowIlPatch[u];

				gsl_vector_scale(x, lhoodFg);


				// Accumulate:
				gsl_vector_add(m, x);


				totalWeight += lhoodFg;
			}
		}

		aRowI += ws;
		aRowImask += ws;
		aRowIlPatch += U;
	}
	gsl_vector_scale(m, 1.0/totalWeight);

	//cout << "UV=" << U*V << "	nPxls=" << nPxls << endl;



	
	//------------------//
	aRowI = (unsigned char*) IpotentialTarget->imageData;
	aRowImask = (unsigned char*) IpotentialTargetMask->imageData;
	aRowIlPatch = (float*) IlPatch->imageData;

	gsl_matrix_set_zero(S);
	//gsl_matrix_set_zero(Sinv);
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			//if((int)aRowImask[u] == 0)
			//if(aRowIlPatch[u]>0.5)
			{

				int uTimes3 = u*3;
				double B = (double) aRowI[uTimes3+0];
				double G = (double) aRowI[uTimes3+1];
				double R = (double) aRowI[uTimes3+2];

			
				eMat->data[0] = B - m->data[0];
				eMat->data[1] = G - m->data[1];
				eMat->data[2] = R - m->data[2];

				gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 , eMat, eMat, 0.0, ee);

				float lhoodFg = 1.0-aRowIlPatch[u];
				gsl_matrix_scale(ee, lhoodFg);
				gsl_matrix_add(S, ee);
			}
		}

		aRowI += ws;
		aRowImask += ws;
		aRowIlPatch += U;

	}

	//gsl_matrix_scale(S, 1.0/(nPxls-1));
	gsl_matrix_scale(S, 1.0/(totalWeight-1));
	//*/


	int s;
	gsl_matrix_memcpy(SLU, S);
	gsl_linalg_LU_decomp(SLU, p3, &s);
	double Sdet = gsl_linalg_LU_det(SLU, s);


		
	int bTarget = 1;
	//if(Sdet<5e6)
	if(Sdet<1e4) // non diag
	{
		bTarget = 0;

		/*
		printGslMatrix(S);
		cout << "Sdet = " << Sdet << endl;
		cvShowImage("Not a Target", IpotentialTarget);
		cvWaitKey();
		*/
	}


	//cvShowImage("Not a Target", IpotentialTarget);


	gsl_vector_free(m);
	gsl_vector_free(x);
	gsl_matrix_free(S);
	gsl_matrix_free(eMat);
	gsl_matrix_free(ee);
	gsl_matrix_free(SLU);
	gsl_permutation_free(p3);




	cvReleaseImage(&IpotentialTarget);
	cvReleaseImage(&IpotentialTargetMask);
	cvReleaseImage(&IsAllegedTarg);
	cvReleaseImage(&IlPatch);





	return bTarget;
}



int drawBoundingBoxes(IplImage* Ilhood, IplImage* ImLandInSea, IplImage* I, IplImage* Ibgr, IplImage* Iboxes, int iFrame, t_track aTrack[])
{
	int U = I->width;
	int V = I->height;
	int lineThickness = 2;
	int flagsForFloodFill = 4| CV_FLOODFILL_MASK_ONLY| CV_FLOODFILL_FIXED_RANGE | (255<<8);
	CvScalar colBlack = cvScalarAll(0);
	CvScalar colWhite = cvScalarAll(255);
	IplImage* Icpy = cvCloneImage(I);
	IplImage* Imask = cvCreateImage(cvSize(U+2,V+2), IPL_DEPTH_8U, 1);
	IplImage* ImaskCorrectSize = cvCreateImage(cvSize(U,V), IPL_DEPTH_8U, 1);
	
	cvCopy(I, Icpy);
	unsigned char* aRowIcpy = (unsigned char*) Icpy->imageData;
	int ws = Icpy->widthStep;

	int nTargets = 0;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRowIcpy[u] == (unsigned char) 0)
			{
				cvZero(Imask);
				cvFloodFill(Icpy, cvPoint(u, v), colWhite, colBlack, colBlack, NULL, flagsForFloodFill, Imask);
				cvSetImageROI(Imask, cvRect(1, 1, U, V));
			
				cvSet(Icpy, colWhite, Imask);
				CvPoint tl, br;
				cvCopy(Imask, ImaskCorrectSize);
				findBoundingBoxDimensions(ImaskCorrectSize, lineThickness, tl, br);

				/*
				cvShowImage("I", I);
				cvShowImage("Imask", Imask);
				cvShowImage("ImaskCorrectSize", ImaskCorrectSize);
				cvShowImage("Icpy", Icpy);
				cvWaitKey();
				//*/

				if(targetAnalysis(Ilhood, ImaskCorrectSize, ImLandInSea, Ibgr, tl, br))
				{

					// Draw the bounding box:
					cvRectangle(Iboxes, tl, br, cvScalar(0,255,255), lineThickness);

					/*
					aTrack[nTargets].aU[iFrame] = ((double)(tl.x + br.x))/2;
					aTrack[nTargets].aV[iFrame] = ((double)(tl.y + br.y))/2;
					*/
					++nTargets;
				}
				else
				{
					//cvRectangle(Iboxes, tl, br, cvScalar(255,100,00), 2);
				}


				cvResetImageROI(Imask);


			} 

		} //u

		aRowIcpy += ws;

	} //v


	//if(nTargets>0)
		//tidyTrackArray(nTargets, iFrame, aTrack);


	cvReleaseImage(&Icpy);
	cvReleaseImage(&Imask);
	cvReleaseImage(&ImaskCorrectSize);

	return nTargets;

}




int postprocess(t_horizon* pHorizon, IplImage* Ilhood, IplImage* ImLandInSea, IplImage* ItoProcess, t_monster* pMons, IplConvKernel* kernel, IplImage* I32F, IplImage* IbgrResults, IplImage* Isignif, int iFrame, t_track aTrack[])
{
	// 0) Usefuls:
	IplImage* I0 = pMons->I1chan8U; 
	IplImage* I1 = pMons->I1chan8U2; 
	IplImage* I2 = pMons->I1chan8U3; 
	IplImage* Iboxes = pMons->I3chan8U ;
	IplImage* IboxesRot = pMons->I3chan8U2 ;

	//cvShowImage("IforSegmo?", I32F);

	cvScale(I32F, I32F, -1.0);
	cvScale(I32F, I32F, 1.0, 1.0);
	cvThreshold(I32F, I1, 0.6, 255, CV_THRESH_BINARY_INV); //.6 //TTT
	cvThreshold(I32F, I0, 0.5, 255, CV_THRESH_BINARY_INV); // .5 
	//cvScale(I32F, I32F, -1.0);
	//cvScale(I32F, I32F, 1.0, 1.0);
#ifdef F_DISPLAY_INTERMEDIATES
	//cvShowImage("4) Thresholded", I1);
	//cvShowImage("ThreGene", I0);
#endif



	cvDilate(I1, I2, kernel);
//#ifdef F_DISPLAY_INTERMEDIATES
	//cvShowImage("Morphoed", I2);
//#endif	
	


	//cvShowImage("5) I0", I0);
	//cvShowImage("5) I2", I2);
	keepSignificantClusters(pMons, I0, I2, Isignif);
//#ifdef F_DISPLAY_INTERMEDIATES
	//cvShowImage("5) Significants", Isignif);
	//displayUnrotatedImage("5) Significants", Isignif);
//#endif




	// 5) Boxes:
	cvSetZero(Iboxes);
	int nTargets = drawBoundingBoxes(Ilhood, ImLandInSea, Isignif, ItoProcess, Iboxes, iFrame, aTrack);
#ifdef F_DISPLAY_INTERMEDIATES
	cvShowImage("Boxes", Iboxes);
#endif

	//cvWarpAffine(Iboxes, IboxesRot, pHorizon->affMat, CV_WARP_INVERSE_MAP|CV_WARP_FILL_OUTLIERS, cvScalarAll(0));
	cvWarpAffine(Iboxes, IboxesRot, pHorizon->affMatRot, CV_WARP_INVERSE_MAP|CV_WARP_FILL_OUTLIERS, cvScalarAll(0));
#ifdef F_DISPLAY_INTERMEDIATES
	cvShowImage("BoxesRot", IboxesRot);
#endif

	cvCopy(IboxesRot, IbgrResults, IboxesRot);
//#ifdef F_DISPLAY_INTERMEDIATES
	cvSetImageROI(IbgrResults, cvRect(pHorizon->xBoxLeft, pHorizon->yBoxTop, pHorizon->U, pHorizon->V));
	cvShowImage("Targets", IbgrResults);
	cvResetImageROI(IbgrResults);
//#endif	



	return nTargets;
}



void checkDeterminants(t_bg* pBg, IplImage* Is)
{
	// Is is silhouette.
	// ImSky is mask, same size.

	cout << "Checking ... " << endl;

	int U = Is->width;
	int V = Is->height;
	unsigned char* aRowIs = (unsigned char*) Is->imageData;
	unsigned char* aRowImSky = (unsigned char*) pBg->ImaskSky->imageData;
	int ws = Is->widthStep;

	int iPxl=0;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRowIs[u]==(unsigned char) 0 && v<0.2*V)
			{
				//cout << pBg->aRGBdiag[iPxl].Sdet << endl;
				printGslMatrix(pBg->aRGBdiag[iPxl].S); 
			}

			++iPxl;
		}

		aRowIs += ws;

	}

	//cvWaitKey();

}




void equaliseImage(CvRect rect, IplImage* Ibgr, IplImage* IbgrEqual)
{

	int U = Ibgr->width;
	int V = Ibgr->height;
	IplImage* I1 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	IplImage* I2 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	IplImage* I3 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	IplImage* I1eq = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	IplImage* I2eq = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	IplImage* I3eq = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	


	cvSplit(Ibgr, I1, I2, I3, NULL);
	cvEqualizeHist(I1, I1eq);
	cvEqualizeHist(I2, I2eq);
	cvEqualizeHist(I3, I3eq);
	//cvMerge(I1eq, I2eq, I3eq, NULL, IbgrEqual);
	cvMerge(I1eq, I2, I3, NULL, IbgrEqual);
	cvShowImage("equalised", IbgrEqual);


	cvCopy(Ibgr, IbgrEqual);

	cvSetImageROI(Ibgr, rect);
	cvSetImageROI(IbgrEqual, rect);
	cvSetImageROI(I1, rect);
	cvSetImageROI(I2, rect);
	cvSetImageROI(I3, rect);
	cvSetImageROI(I1eq, rect);
	cvSetImageROI(I2eq, rect);
	cvSetImageROI(I3eq, rect);
	cvSplit(Ibgr, I1, I2, I3, NULL);
	cvEqualizeHist(I1, I1eq);
	cvEqualizeHist(I2, I2eq);
	cvEqualizeHist(I3, I3eq);
	//cvMerge(I1eq, I2eq, I3eq, NULL, IbgrEqual);	
	cvMerge(I1eq, I2, I3, NULL, IbgrEqual);	
	cvShowImage("I1", I1);
	cvShowImage("I1eq", I1eq);
	cvShowImage("I2", I2);
	cvShowImage("I2eq", I2eq);
	cvShowImage("I3", I3);
	cvShowImage("I3eq", I3eq);

	cvResetImageROI(Ibgr);
	cvResetImageROI(IbgrEqual);


	cvReleaseImage(&I1);
	cvReleaseImage(&I2);
	cvReleaseImage(&I3);
	cvReleaseImage(&I1eq);
	cvReleaseImage(&I2eq);
	cvReleaseImage(&I3eq);

	
}




void calculateImageStatistics(IplImage* I)
{
	// 8-bit, 3-channel.

	int U = I->width;	
	int V = I->height;	
	unsigned char* aRow = (unsigned char*) I->imageData;



	double cumulB = 0;
	double cumulBB = 0;
	double cumulG = 0;
	double cumulGG = 0;
	double cumulR = 0;
	double cumulRR = 0;

	for(int iCh=0; iCh<1; ++iCh)
	{

		for(int v=0; v<V; ++v)
		{
			for(int u=0; u<U; ++u)
			{
				// Scanning for black pixels:
				int uTimes3 = u*3;
				double B = (double) aRow[uTimes3+0] ;
				double G = (double) aRow[uTimes3+1] ; 
				double R = (double) aRow[uTimes3+2] ; 

				cumulB	+= B;
				cumulBB += B*B;
				cumulG	+= G;
				cumulGG += G*G;
				cumulR	+= R;
				cumulRR += R*R;

			}
		
			aRow += I->widthStep;

		}

	}

	int N = U*V;
	double meanB = cumulB/N;
	double varB = cumulBB/N - meanB*meanB;
	double stdB = sqrt(varB);

	double meanG = cumulG/N;
	double varG = cumulGG/N - meanG*meanG;
	double stdG = sqrt(varG);
	
	double meanR = cumulR/N;
	double varR = cumulRR/N - meanR*meanR;
	double stdR = sqrt(varR);

	double cumul = cumulB + cumulG + cumulR;
	double cumulSq = cumulBB + cumulGG + cumulRR;
	double mean = cumul/(3*N);
	double var = cumulSq/(3*N) - mean*mean;
	double std = sqrt(var);

	/*
	cout << "mean = (" << meanB << ", " << meanG << ", " << meanR << ")" << endl;
	cout << "std = (" << stdB << ", " << stdG << ", " << stdR << ")" << endl;
	*/
	cout << "mean = " << mean << "	std = " << std << endl;
}




void verifyImage(IplImage* I)
{
	// I is 8-bit 3-channel.

	int U = I->width;	
	int V = I->height;	
	unsigned char* aRow = (unsigned char*) I->imageData;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			int uTimes3 = u*3;
			unsigned char& chB = aRow[uTimes3+0];
			unsigned char& chG = aRow[uTimes3+1];
			unsigned char& chR = aRow[uTimes3+2];
			



		
			if((int) chB > 255)
				chB = (unsigned char) 255;
			if((int) chB < 0)
				chB = (unsigned char) 0;

			if((int) chG > 255)
				chG = (unsigned char) 255;
			if((int) chG < 0)
				chG = (unsigned char) 0;

			if((int) chR > 255)
				chR = (unsigned char) 255;
			if((int) chR < 0)
				chR = (unsigned char) 0;
		
		}

		aRow += I->widthStep;
	}
}


void linearTransform(IplImage* I, double alpha, double beta)
{
	// I is 8bit 3chan.

	
	
	CvSize size = cvSize(I->width, I->height);
	IplImage* Isc = cvCreateImage(size, IPL_DEPTH_32F, 3);
	double mean = 160.0/255;
	cvConvertScale(I, Isc, 1.0/255); 
	cvConvertScale(Isc, Isc, 1.0, -mean);
	cvConvertScale(Isc, Isc, 1.0, beta+mean); // alpha
	cvConvertScale(Isc, I, 255.0, 0);
	verifyImage(I);
	cvShowImage("Scaled", I);
	cvReleaseImage(&Isc);
	
}



int isWhite(IplImage* I)
{
	int U = I->width;	
	int V = I->height;	
	unsigned char* aRow = (unsigned char*) I->imageData;
	int ws = I->widthStep;

	int bWhite = 1;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRow[u] != (unsigned char) 255)
				bWhite = 0;
		}

		aRow += ws;
	}

	return bWhite;
}



void finalAdaptations(IplImage* IsSignif, IplImage* Ilhood, IplImage* IlImproved2)
{

	
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);
	cvErode(IsSignif, IsSignif);

	int U = IsSignif->width;	
	int V = IsSignif->height;	
	unsigned char* aRowSign = (unsigned char*) IsSignif->imageData;
	float* aRowLhood = (float*) Ilhood->imageData;
	float* aRowLImpr = (float*) IlImproved2->imageData;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			aRowLImpr[u] = 1- (1-aRowLhood[u]) * (1-((float) ((int) aRowSign[u]))/255.0 );

		}

		aRowSign += IsSignif->widthStep;
		aRowLhood += Ilhood->widthStep;
		aRowLImpr += IlImproved2->widthStep;
	}


	if(isWhite(IsSignif))
		;
	else
		cvCopy(IlImproved2, Ilhood);

	cvShowImage("Improved2", Ilhood);



}


