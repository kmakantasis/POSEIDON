// Everything to do with smart convo functions and speedups.




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
#include "foregroundconstants.h"

#include <opencv2\video\background_segm.hpp>


using std::cout;
using std::endl;

extern t_settings* g_pSetts;



float multiply(float a, float b)
{
	return a*b;
}



float pcvMultiplySum(int U, int u0, int v0, int K, int k, int kEff, IplImage* I1, IplImage* I2)
{
	// Idst = I1*I2, pixelwise.


	// 1) Usefuls:
	int kMinuskEff = k-kEff;
	int kPluskEff = k+kEff;
	float* aRowI2 =	(float*) I2->imageData + kMinuskEff*K;
	float* aRowI1 = (float*) I1->imageData + (v0+kMinuskEff)*U + u0;

	
	float sumVal = 0.0;
	for(int v=kMinuskEff; v<=kPluskEff; ++v)
	{
		for(int u=kMinuskEff; u<=kPluskEff; ++u)
		{
			//sumVal += multiply(aRowI1[u], aRowI2[u]);
			sumVal += aRowI1[u] * aRowI2[u];
		}

		aRowI1 += U ;
		aRowI2 += K ;
	}

	return sumVal;
}




float pcvMultiplySumRect(int u0, int v0, int KeffU, int KeffV, IplImage* I, IplImage* Ik)
{
	// Idst = I1*I2, pixelwise.

	// 1) Usefuls:
	int U = I->width;
	int KU = Ik->width;
	int kU = (KU-1)/2;
	int kEffU = (KeffU-1)/2;	

	int KV = Ik->height;
	int kV = (KV-1)/2;
	int kEffV = (KeffV-1)/2;	
	
	float* __restrict aRowI =  (float*)  I->imageData + (v0+kV-kEffV)*U + u0;
	float* __restrict aRowIk =	(float*) Ik->imageData + (kV-kEffV)*KU;

	
	float sumVal = 0.0;
	for(int v=kV-kEffV; v<=kV+kEffV; ++v)
	{
		for(int u=kU-kEffU; u<=kU+kEffU; ++u)
		{
			//sumVal += multiply(aRowI[u], aRowIk[u]);
			sumVal += aRowI[u] * aRowIk[u];
		}

		aRowI += U ;
		aRowIk += KU ;
	}

	return sumVal;
}





float pTotalKernelIntensityRect(int KeffU, int KeffV, IplImage* Ik)
{

	int kEffU = (KeffU-1)/2;	
	int KU = Ik->width;
	int kU = (KU-1)/2;

	int kEffV = (KeffV-1)/2;	
	int KV = Ik->height;
	int kV = (KV-1)/2;


	float sumPix=0.0;
	for(int v=kV-kEffV; v<=kV+kEffV; ++v)
	{
		float* aRowI = (float*) Ik->imageData + KU*v;
		
		for(int u=kU-kEffU; u<=kU+kEffU; ++u)
			sumPix += aRowI[u];
	}

	return sumPix;
}




void pRecalculateKernelRect(int KeffU, int KeffV, double sigma, IplImage* Ik)
{	
	// For kernel images, 32bit, 1chan, KU by KV, each ODD.

	int kEffU = (KeffU-1)/2;	
	int KU = Ik->width;
	int kU = (KU-1)/2;

	int kEffV = (KeffV-1)/2;	
	int KV = Ik->height;
	int kV = (KV-1)/2;

	double sigmaSqTimes2 = 2*sigma*sigma;

	for(int v=kV-kEffV; v<=kV+kEffV; ++v)
	{
		float* aRowIk = (float*) Ik->imageData + KU*v;
	
		double errV = v-kV;
		double errVsq = errV*errV;

		for(int u=kU-kEffU; u<=kU+kEffU; ++u)
		{
			double errU = u-kU;
			double mahalanoSq = (errU*errU + errVsq)/sigmaSqTimes2;
			aRowIk[u] = exp(-mahalanoSq);
		}

	}
	
	//cvShowImage("kernel", Ik);
	//cvWaitKey();
	float totalIntensity = pTotalKernelIntensityRect(KeffU, KeffV, Ik);
	cvScale(Ik, Ik, 1/totalIntensity);

}







void pSmartSmoothNewFast(t_monster* pMons, double yMh, IplImage* ImAvoid, IplImage* I, IplImage* Iout)
{
	// Convolve with Gaussian whose standard deviation is sigma, 
	// where sigma is a function of distance from a given line.

	int U = I->width;
	int V = I->height;
	int Vunrot = g_pSetts->V;
	double scFctr = ((double) Vunrot)/180.0;
	double lamda = g_pSetts->DeltaStDev/sqrt((double)Vunrot) * scFctr;
	double sigmaMin = g_pSetts->sigmaMin;
	double thSigma = g_pSetts->thSigma;

	// Kernel stuff:
	int K = pMons->K;
	int k = pMons->k;



	int vHorizon = (int) (yMh+0.5);
	


	
	cvSet(Iout, cvScalarAll(1));
	
	// Loop:
	for(int iRun=0; iRun<2; ++iRun)
	{

		unsigned char* aRowImAvoid = (unsigned char*) ImAvoid->imageData;
		double sigma = sigmaMin;
		double sigmaPrev=0;
		int KeffV=1, KeffU=1, kEffV=0, kEffU=0;
		int& Keff = (iRun==0)? KeffU : KeffV;
		int& kEff = (iRun==0)? kEffU : kEffV;
		IplImage* Ik = pMons->aIk[iRun];
		float* aRowIout = (float*) Iout->imageData;
		IplImage* Iin = pMons->aIin[iRun];
		CvPoint offset = (iRun==0)? cvPoint(k,0) : cvPoint(0,k);
		int Uin = Iin->width;
		IplImage* IcopyFrom = (iRun==0)? I: Iout;
		cvCopyMakeBorder(	IcopyFrom, 
							Iin, 
							offset, 
							IPL_BORDER_CONSTANT, 
							cvScalarAll(1)	// white border
						);

		for(int v=0; v<V; ++v)
		{

			// Calculate sigma:
			double d = ((float) (v - vHorizon));

			if(d>0)
				sigma = lamda*sqrt(d) + sigmaMin;
				
		
			if(fabs(sigma-sigmaPrev) > thSigma)
			{
				kEff = (int) sqrt(-2*sigma*sigma*log(thSigma) );
				Keff = 2*kEff+1;


				pRecalculateKernelRect(KeffU, KeffV, sigma, Ik);

				sigmaPrev = sigma;

				//cout << sigma << endl;
			}



			for(int u=0; u<U; ++u)
			{
				if(aRowImAvoid[u] == (unsigned char) 0)
					aRowIout[u] = pcvMultiplySumRect( u, v, KeffU, KeffV, Iin, Ik);
			}

			aRowIout += U;
			aRowImAvoid += ImAvoid->widthStep;

		} // for v

	} // for iRun


}




