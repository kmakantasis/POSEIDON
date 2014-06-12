// background model class implementation
//
// Paris Kaimakis 19/10/12




#include <iostream>
#include <cfloat>
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

#include <opencv2/video/background_segm.hpp>


using std::cout;
using std::endl;

extern t_settings* g_pSetts;


// One of the following must be on:
//#define LF_FULL_COVAR_TOO
#define LF_FULL_COVAR_TOO_1


//#define F_DISPLAY_INTERMEDIATES


void createBg(t_bg* pBg, int U, int V)
{
	cout << "\tCreating background class... ";
	cout.flush();

	pBg->U = U;
	pBg->V = V;


	int nGaussians = U*V;
#ifdef LF_FULL_COVAR_TOO
	pBg->aRGB = (t_nrml*) calloc(nGaussians, sizeof(t_nrml));
#endif
#ifdef LF_FULL_COVAR_TOO_1
	pBg->aRGBdiag = (t_nrml*) calloc(nGaussians, sizeof(t_nrml));
#endif
	for(int i=0; i<nGaussians; ++i)
	{
#ifdef LF_FULL_COVAR_TOO
		createNormal(&pBg->aRGB[i], 3);
#endif
#ifdef LF_FULL_COVAR_TOO_1
		createNormal(&pBg->aRGBdiag[i], 3);
#endif
	}


	CvSize size = cvSize(U, V);
	for(int i=0; i<NUM_FRAMES_FOR_TRAINING; ++i)
		pBg->aIsample[i] = cvCreateImage(size, IPL_DEPTH_8U, 3);


#ifdef LF_FULL_COVAR_TOO
	pBg->Ilhood = cvCreateImage(size, IPL_DEPTH_32F, 1);
#endif
#ifdef LF_FULL_COVAR_TOO_1
	pBg->IlhoodDiag = cvCreateImage(size, IPL_DEPTH_32F, 1);
#endif

	pBg->Imean = cvCreateImage(size, IPL_DEPTH_32F, 3);
	pBg->Isdet = cvCreateImage(size, IPL_DEPTH_32F, 1);
	pBg->ImaskSdetZero = cvCreateImage(size, IPL_DEPTH_8U, 1);
	pBg->ImaskRubbish = cvCreateImage(size, IPL_DEPTH_8U, 1);
	
	cvSetZero(pBg->ImaskSdetZero );
	cvSet(pBg->ImaskRubbish, cvScalarAll(255));
	cvRectangle(pBg->ImaskRubbish, cvPoint(0,0), cvPoint(U/5, V/5), cvScalarAll(0), -1);
	cvRectangle(pBg->ImaskRubbish, cvPoint(U/5, 9*V/10), cvPoint(4*U/5, V-1), cvScalarAll(0), -1);
	cvRectangle(pBg->ImaskRubbish, cvPoint(0, 0), cvPoint(U/2, V/10), cvScalarAll(0), -1);

	// Monster:
	createForegroundMonster(U, V, &pBg->mons);




	pBg->mog.set("history", 40);
	pBg->mog.set("nmixtures", 5);
	pBg->mog.set("backgroundRatio", 0.25);
	pBg->mog.set("noiseSigma", .01);



	pBg->ImaskSky = cvCreateImage(size, IPL_DEPTH_8U, 1);
	cvSetZero(pBg->ImaskSky);

	int Vunrot = g_pSetts->V;
	int kerSize = (int) (3.0);// /180*Vunrot); //3.0
	int kerAnch = kerSize/2;
	pBg->kernel = cvCreateStructuringElementEx(kerSize, kerSize, kerAnch, kerAnch, CV_SHAPE_RECT);

	pBg->IforSegmo = cvCreateImage(size, IPL_DEPTH_32F, 1);



	pBg->ImaskSea = cvCreateImage(size, IPL_DEPTH_8U, 1);
	cvSetZero(pBg->ImaskSea);
	cvCopy(pBg->ImaskSea, pBg->mons.ImSea);

	pBg->IsSignificant = cvCreateImage(size, IPL_DEPTH_8U, 1);
	pBg->IlImproved2 = cvCreateImage(size, IPL_DEPTH_32F, 1);



	pBg->Iresized = cvCreateImage(size, IPL_DEPTH_8U, 3);
	pBg->Iresults = cvCreateImage(size, IPL_DEPTH_8U, 3);
	pBg->ItoProcess = cvCreateImage(size, IPL_DEPTH_8U, 3);
	
	pBg->ImLandInSea = cvCreateImage(size, IPL_DEPTH_8U, 1);

	cvSetZero(pBg->Iresized);
	cvSetZero(pBg->Iresults);
	cvSetZero(pBg->ItoProcess);



	createStuffForPaper(&pBg->paper);


	// Output to be used to form the posterior:
	pBg->IjointBG = cvCreateImage(cvSize(g_pSetts->U, g_pSetts->V), IPL_DEPTH_32F, 1);


	cout << "Done!" << endl;


}





void destroyBg(t_bg* pBg)
{
	
	cout << "\tDestroying background class... " ;
	cout.flush();

	
	// Gaussians:
	int nGaussians = pBg->U*pBg->V;
	for(int i=0; i<nGaussians; ++i)
	{
#ifdef LF_FULL_COVAR_TOO
		deleteNormal(&pBg->aRGB[i]);
#endif
#ifdef LF_FULL_COVAR_TOO_1
		deleteNormal(&pBg->aRGBdiag[i]);
#endif
	}
#ifdef LF_FULL_COVAR_TOO
	free(pBg->aRGB);
#endif
#ifdef LF_FULL_COVAR_TOO_1
	free(pBg->aRGBdiag);
#endif

	for(int i=0; i<NUM_FRAMES_FOR_TRAINING; ++i)
		cvReleaseImage(&pBg->aIsample[i]);
	

#ifdef LF_FULL_COVAR_TOO
	cvReleaseImage(&pBg->Ilhood);
#endif
#ifdef LF_FULL_COVAR_TOO_1
	cvReleaseImage(&pBg->IlhoodDiag);
#endif

	cvReleaseImage(&pBg->Imean);
	cvReleaseImage(&pBg->Isdet);
	cvReleaseImage(&pBg->ImaskSdetZero);
	cvReleaseImage(&pBg->ImaskRubbish);
	cvReleaseImage(&pBg->ImaskSea);
	cvReleaseImage(&pBg->ImaskSky);
	cvReleaseImage(&pBg->IforSegmo);
	cvReleaseImage(&pBg->IsSignificant);
	cvReleaseImage(&pBg->IlImproved2);
	cvReleaseImage(&pBg->ImLandInSea);


	destroyForegroundMonster(&pBg->mons);


	cvReleaseStructuringElement(&pBg->kernel);


	cvReleaseImage(&pBg->Iresized);
	cvReleaseImage(&pBg->Iresults);
	cvReleaseImage(&pBg->ItoProcess);


	deleteStuffForPaper(&pBg->paper);

	// Output to be used to form the posterior:
	cvReleaseImage(&pBg->IjointBG);

	cout << "Done!" << endl;

}







void trainBgModel(CvCapture* capture, IplImage* Iframe, t_horizon* pHorizon, IplImage* Iresized, CvMat* affMat, IplImage* ItoProcess, t_bg* pBg)
{
	cout << "\tTraining background model... " ;
	cout.flush();


	recordBgTrainingSet(capture, Iframe, pHorizon, Iresized, affMat, ItoProcess, pBg);
#ifdef LF_FULL_COVAR_TOO
	calculateBgModel(F_NONDIAG, pBg);
#endif
#ifdef LF_FULL_COVAR_TOO_1
	calculateBgModel(F_DIAG, pBg);
#endif

	cout << "Done!" << endl;

}


void recordBgTrainingSet(CvCapture* capture, IplImage* Iframe, t_horizon* pHorizon, IplImage* Iresized, CvMat* affMat, IplImage* ItoProcess, t_bg* pBg)
{


	for(int i=0; i<NUM_FRAMES_FOR_TRAINING; ++i)
	{
		
		Iframe = cvQueryFrame(capture);
		cvResize(Iframe, Iresized);

		makeHorizonHorizontalPaddedWay(pHorizon, Iresized, pBg->Iresized, pBg->ItoProcess);
		cvCopy(pBg->ItoProcess, pBg->Iresized);
		cvCopy(pBg->ItoProcess, pBg->Iresults);



		//cvCopyImage(ItoProcess, pBg->aIsample[i]);
		cvCopy(ItoProcess, pBg->aIsample[i]);

		/*
#ifdef F_PAPER
		displayPixelIntensityGraphsForPaper(Iresized, i, &pBg->paper);
#endif
		*/
	}


	//cvWaitKey();
}


void trainPixelCovarsNew(int flag, t_bg* pBg)
{
	int U = pBg->U;
	int V = pBg->V;
	int nPxls = U*V;



	for(int iPxl=0; iPxl<nPxls; ++iPxl)
	{
		t_nrml* aGauss = pBg->aRGBdiag;

		gsl_vector* m = aGauss[iPxl].m;
		gsl_matrix* S = aGauss[iPxl].S;
		double& Sdet = aGauss[iPxl].Sdet;
		gsl_matrix* Sinv = aGauss[iPxl].Sinv;

		double var = 100;
		gsl_matrix_set_identity(S);
		gsl_matrix_scale(S, var);
		gsl_matrix_set_identity(Sinv);
		gsl_matrix_scale(Sinv, 1/var);
		Sdet = var*var*var;

	}
}



void calculateBgModel(int flag, t_bg* pBg)
{
	trainPixelMeans(flag, pBg);
	trainPixelCovars(flag, pBg);
}






void trainPixelMeans(int flag, t_bg* pBg)
{

	int U = pBg->U;
	int V = pBg->V;
	int nPxls = U*V;
	gsl_vector* x = gsl_vector_calloc(3);
	IplImage* Imask = pBg->ImaskRubbish;

	t_nrml* aGauss = NULL;
	switch(flag)
	{
	case F_NONDIAG:
		aGauss = pBg->aRGB;
		break;
	case F_DIAG:
		aGauss = pBg->aRGBdiag;
		break;
	}



	for(int iPxl=0; iPxl<nPxls; ++iPxl)
	{
		gsl_vector_set_zero(aGauss[iPxl].m);
	}




	int nFrames = 0;
	for(int ifr=0; ifr<NUM_FRAMES_FOR_TRAINING; ++ifr)
	{
		IplImage* Is = pBg->aIsample[ifr];
	
		int iPxl = 0;
		for(int v=0; v<V; ++v)
		{
			unsigned char* aRowI = (unsigned char*) Is->imageData + Is->widthStep*v;
			unsigned char* aRowImask = (unsigned char*) Imask->imageData + Imask->widthStep*v;



			for(int u=0; u<U; ++u)
			{
				//if(aRowImask[u] == (unsigned char)255)
				{
					int uTimes3 = u*3;
					double B = (double) aRowI[uTimes3+0] ;
					double G = (double) aRowI[uTimes3+1] ; 
					double R = (double) aRowI[uTimes3+2] ; 

					x->data[0] = B;
					x->data[1] = G;
					x->data[2] = R;


					gsl_vector* m = aGauss[iPxl].m;
					gsl_vector_add(m, x);
				}
				++iPxl;
			} //u
		} //v

		++nFrames;
	} //ifr



	for(int iPxl=0; iPxl<nPxls; ++iPxl)
		gsl_vector_scale(aGauss[iPxl].m, 1.0/nFrames);



}




void trainPixelCovars(int flag, t_bg* pBg)
{
	int U = pBg->U;
	int V = pBg->V;
	int nPxls = U*V;
	gsl_matrix* eMat = pBg->mons.eMat;
	gsl_matrix* ee = pBg->mons.ee;
	gsl_matrix* SLU = pBg->mons.SLU;
	gsl_permutation* p3 = pBg->mons.p3;

	t_nrml* aGauss = NULL;
	switch(flag)
	{
	case F_NONDIAG:
		aGauss = pBg->aRGB;
		break;
	case F_DIAG:
		aGauss = pBg->aRGBdiag;
		break;
	}





	for(int iPxl=0; iPxl<nPxls; ++iPxl)
		gsl_matrix_set_zero(aGauss[iPxl].S);




	int nFrames = 0;
	for(int ifr=0; ifr<NUM_FRAMES_FOR_TRAINING; ++ifr)
	{
		IplImage* Is = pBg->aIsample[ifr];
	
		int iPxl = 0;
		for(int v=0; v<V; ++v)
		{
			unsigned char* aRowI = (unsigned char*) Is->imageData + Is->widthStep*v;

			for(int u=0; u<U; ++u)
			{
				// Usefuls:
				gsl_vector* m = aGauss[iPxl].m;
				gsl_matrix* S = aGauss[iPxl].S;


				int uTimes3 = u*3;
				double B = (double) aRowI[uTimes3+0] ;
				double G = (double) aRowI[uTimes3+1] ; 
				double R = (double) aRowI[uTimes3+2] ; 

				eMat->data[0] = B - m->data[0];
				eMat->data[1] = G - m->data[1];
				eMat->data[2] = R - m->data[2];

				gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, eMat, eMat, 0.0, ee);
				
				gsl_matrix_add(S, ee);


				++iPxl;
			} //u
		} //v

		++nFrames;
	} //ifr




	cvZero(pBg->ImaskSdetZero);
	for(int iPxl=0; iPxl<nPxls; ++iPxl)
	{
		gsl_matrix* S = aGauss[iPxl].S;
		double& Sdet = aGauss[iPxl].Sdet;
		gsl_matrix* Sinv = aGauss[iPxl].Sinv;

		gsl_matrix_scale(S, 1.0/((double)(nFrames-1)));
		
		if(flag == F_DIAG)
		{
			S->data[1] = 0;
			S->data[2] = 0;
			S->data[3] = 0;
			S->data[5] = 0;
			S->data[6] = 0;
			S->data[7] = 0;
		}


		int s;
		gsl_matrix_memcpy(SLU, S);
		gsl_linalg_LU_decomp(SLU, p3, &s);
		Sdet = gsl_linalg_LU_det(SLU, s);




		


		removeSaturationInducedSingularities(&pBg->mons, aGauss[iPxl].m, 0.111, aGauss[iPxl].Sdet, aGauss[iPxl].S);

		
		Sdet = aGauss[iPxl].Sdet;

		
		double tSdet = 1;//100;// SMALL_VALUE; //1; //FLT_MIN;
		double f = 1; // 1*std::exp(-Sdet/tSdet);
		//cout << "f = " << f << endl;
		double fMagn = f*1.2;
		double fAttn = 0*0.8;
		double fMin = 1;
		//double 
		int iC = 0;

		//while(0)
		while(Sdet<tSdet)
		{
			
			double& s00 = S->data[0];
			double& s11 = S->data[4];
			double& s22 = S->data[8];



			S->data[1] *= fAttn;
			S->data[2] *= fAttn;
			S->data[3] *= fAttn;
			S->data[5] *= fAttn;
			S->data[6] *= fAttn;
			S->data[7] *= fAttn;

			S->data[0] *= fMagn;
			S->data[4] *= fMagn;
			S->data[8] *= fMagn;
			
			S->data[0] = std::max(S->data[0], fMin);
			S->data[4] = std::max(S->data[4], fMin);
			S->data[8] = std::max(S->data[8], fMin);


			// Recalculate Sdet:
			gsl_matrix_memcpy(SLU, S);
			gsl_linalg_LU_decomp(SLU, p3, &s);
			Sdet = gsl_linalg_LU_det(SLU, s);

			++iC;
			
		} //while Sdet < tSdet


		
		// Sinv:
		gsl_linalg_LU_invert(SLU, p3, Sinv);



	}



}



float pixelColourLikelihood(int flag, double B, double G, double R, t_nrml* gaussian)
{
	// Usefuls:
	gsl_matrix* Sinv = gaussian->Sinv;
	gsl_vector* e = gaussian->e;
	gsl_vector* m = gaussian->m;
	gsl_vector* SinvTimese = gaussian->SinvTimese;
	double Sdet = gaussian->Sdet;


	double* eData = e->data;
	double* mData = m->data;
	eData[0] = mData[0] - B;
	eData[1] = mData[1] - G;
	eData[2] = mData[2] - R;


	double mahalanobisSq;
	gsl_blas_dgemv(CblasNoTrans, 1.0, Sinv, e, 0.0, SinvTimese);
	gsl_blas_ddot(e, SinvTimese, &mahalanobisSq);


	return std::exp(-0.5*mahalanobisSq);
}




float pixelColourLikelihoodFast(double B, double G, double R, t_nrml* gaussian)
{
	// Usefuls:
	double* mData = gaussian->m->data;
	double* SinvData = gaussian->Sinv->data;


	// 1) Load e:
	double e1 = mData[0] - B;
	double e2 = mData[1] - G;
	double e3 = mData[2] - R;

	double s11 = SinvData[0];
	double s12 = SinvData[1];
	double s13 = SinvData[2];
	double s22 = SinvData[4];
	double s23 = SinvData[5];
	double s33 = SinvData[8];



	double mahalanoSq = (	+ e1*e1*s11 + 2*e1*e2*s12 + 2*e1*e3*s13 
										+ e2*e2*s22 + 2*e2*e3*s23 
										+ e3*e3*s33 
									 );

	double likelihood = exp(-0.5*mahalanoSq);
	

	/*
	//cout << FLT_MIN << endl;
	cout << likelihood << endl;

	
	if(likelihood>0 && likelihood<FLT_MIN)
		cout << "AAA" << endl;
	//*/

	return likelihood;



	/*
	return std::exp(-0.5*(	+ e1*e1*s11 + 2*e1*e2*s12 + 2*e1*e3*s13 
										+ e2*e2*s22 + 2*e2*e3*s23 
										+ e3*e3*s33 
									 )
							   );	
							   */


}



float pixelColourLikelihoodFastDiag(double B, double G, double R, t_nrml* gaussian)
{
	// Like pixelColourLikelihoodFast but assuming diagonal covariance.

	// Usefuls:
	double* __restrict mData = gaussian->m->data;
	double* __restrict SinvData = gaussian->Sinv->data;

	// 1) Load e:
	double e1 = mData[0] - B;
	double e2 = mData[1] - G;
	double e3 = mData[2] - R;

	double s11 = SinvData[0];
	double s22 = SinvData[4];
	double s33 = SinvData[8];

	return std::exp( -0.5*(	+ e1*e1*s11 
							+ e2*e2*s22 
							+ e3*e3*s33 
						   ) 
				   );
}






void backgroundLikelihood(t_nrml aGaussian[], IplImage* ImAvoid, IplImage* Iframe, IplImage* Ilhood, IplImage* ImaskSky)
{
	// Usefuls:
	int U = Iframe->width;
	int V = Iframe->height;
	int IframeWidthStep = Iframe->widthStep;

	unsigned char* __restrict aRowI = (unsigned char*) Iframe->imageData;
	unsigned char* __restrict aRowImAvoid = (unsigned char*) ImAvoid->imageData;
	float* __restrict aRowIlhood = (float*) Ilhood->imageData;


	cvSet(Ilhood, cvScalarAll(1));

	int iPxl = 0;
	for(int v=0; v<V; ++v)
	{      

		for(int u=0; u<U; ++u)
		{			

			if(aRowImAvoid[u] == (unsigned char) 0)
			{
				// Usefuls:
				t_nrml* __restrict gaussian = &aGaussian[iPxl];

				int uTimes3 = u*3;
				double B = (double) aRowI[uTimes3++];
				double G = (double) aRowI[uTimes3++];
				double R = (double) aRowI[uTimes3];


				// Calculate likelihood for this colour:
#ifdef LF_FULL_COVAR_TOO
				aRowIlhood[u] = pixelColourLikelihoodFast(B, G, R, gaussian);
#endif
#ifdef LF_FULL_COVAR_TOO_1
				aRowIlhood[u] = pixelColourLikelihoodFastDiag(B, G, R, gaussian);
#endif

			} // if avoid-mask

			++iPxl;
		}

		
		aRowImAvoid += ImAvoid->widthStep;
		aRowI += IframeWidthStep;
		aRowIlhood += U;
		
		
	}


#ifdef F_DISPLAY_INTERMEDIATES
	cvShowImage("Likelihood", Ilhood);
#endif

}




void pThresholdToMin(IplImage* Isrc, IplImage* Idst, float thresh)
{

	int V = Isrc->height;
	int U = Isrc->width;
	float* aRowIsrc = (float*) Isrc->imageData;
	float* aRowIdst = (float*) Idst->imageData ;


	for(int v=0; v<V; ++v)
	{

		for(int u=0; u<U; ++u)
		{			
			aRowIdst[u] = (aRowIsrc[u] < thresh) ? thresh : aRowIsrc[u];
		}

		aRowIsrc += U;
		aRowIdst += U;
	}


}



float pTotalImageIntensity(IplImage* I)
{
	// Add all the pixels in an image.
	// 1 channel, 32bit.

	int V = I->height;
	int U = I->width;
	float* aRowI = (float*) I->imageData;

	float sumPix=0.0;
	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
			sumPix += aRowI[u];

		aRowI += U;
	}

	return sumPix;
}



float pTotalKernelIntensity(int Keff, IplImage* I)
{
	// Add all the pixels in an image.
	// 1 channel, 32bit.
	// Kernel images, ie square, odd-sized side.
	// This function may be called while being on a specifi ROI. Hence the Keff, K.
	// 
	// K is the original image size.

	int kEff = (Keff-1)/2;
	int K = I->width;
	int k = (K-1)/2;

	float sumPix=0.0;
	for(int v=k-kEff; v<=k+kEff; ++v)
	{
		float* aRowI = (float*) I->imageData + K*v;
		
		for(int u=k-kEff; u<=k+kEff; ++u)
			sumPix += aRowI[u];
	}

	return sumPix;
}




void pRecalculateKernel(int Keff, double sigma, IplImage* Ik)
{	
	// For kernel images, 32bit, 1chan, square, each side is ODD number.
	//

	int kEff = (Keff-1)/2;	
	int K = Ik->width;
	int k = (K-1)/2;
	double sigmaSqTimes2 = 2*sigma*sigma;

	for(int v=k-kEff; v<=k+kEff; ++v)
	{
		float* aRowIk = (float*) Ik->imageData + K*v;
	
		double errV = v-k;
		double errVsq = errV*errV;

		for(int u=k-kEff; u<=k+kEff; ++u)
		{
			double errU = u-k;
			double mahalanoSq = (errU*errU + errVsq)/sigmaSqTimes2;
			aRowIk[u] = exp(-mahalanoSq);
		}

	}
	
	// Normalise:
	float totalIntensity = pTotalKernelIntensity(Keff, Ik);
	cvScale(Ik, Ik, 1/totalIntensity);

}





void improveBackgroundLikelihood(t_monster* pMons, double yMh, IplImage* ImAvoid, IplImage* IlBack32F, IplImage* IforSegmo, IplImage* IjointBG)
{

	int U = IlBack32F->width;
	int V = IlBack32F->height;
	double scFct = ((double)V)/180;
	double sigma = pMons->sigma;
	int kerSize = pMons->kerSize;
	IplImage* IlBack32Fconvo = pMons->I32F1;
	IplImage* IwithBorder = pMons->IwithBorder;



	
	// Smart-convo according to knowledge of horizon:	
	pSmartSmoothNewFast(pMons, yMh, ImAvoid, IlBack32F, IlBack32Fconvo);	
	//displayUnrotatedImage("4) Horizon-Dependent Convo I", IlBack32Fconvo);

	cvMul(IlBack32F, IlBack32Fconvo, IlBack32F);
	cvMul(IlBack32Fconvo, IlBack32Fconvo, IforSegmo);



	returnUnrotatedImage("Unrotated?", IforSegmo, IjointBG);


	pThresholdToMin(IlBack32F, IlBack32F, 0.01);


#ifdef F_DISPLAY_INTERMEDIATES
	cvShowImage("Improved", IlBack32F);
	cvShowImage("IforSegmo", IforSegmo);
#endif



	//cvShowImage("AAA", IforSegmo);

}




void removeBackground(t_monster* pMons, IplImage* ImAvoid, double yMh, t_nrml aGaussian[], IplImage* ItoProcess, IplImage* Ilhood, IplImage* ImaskSky, IplImage* IforSegmo, IplImage* IjointBG)
{

	

	// 1) Likelihood according to bg-model:
	// Note: Land-in-sea likelihoods are 1.0!:
	backgroundLikelihood(aGaussian, ImAvoid, ItoProcess, Ilhood, ImaskSky);

	// 2) Improve likelihood:
	improveBackgroundLikelihood(pMons, yMh, ImAvoid, Ilhood, IforSegmo, IjointBG);
	//cvShowImage("4) Improved Likelihood", Ilhood);
	//displayUnrotatedImage("4) Horizon-Dependent Convo II", Ilhood);
	//*/




}






void adaptBackground(t_monster* pMons, double alpha, IplImage* I, t_bg* pBg)
{
#ifdef LF_FULL_COVAR_TOO
	adaptBackgroundThisModel(pBg->ImLandInSea, pMons, alpha, I, pBg->Ilhood, pBg->aRGB);
#endif
#ifdef LF_FULL_COVAR_TOO_1	
	adaptBackgroundThisModel(pBg->ImLandInSea, pMons, alpha, I, pBg->IlhoodDiag, pBg->aRGBdiag);
#endif
}



void pGslMatrixMemcpyDiag(gsl_matrix* dst, gsl_matrix* src)
{
	double* __restrict dstData = dst->data;
	double* __restrict srcData = src->data;

	dstData[0] = srcData[0];
	dstData[4] = srcData[4];
	dstData[8] = srcData[8];

}


void pGslMatrixScaleDiag(gsl_matrix* m, double s)
{
	double* __restrict mData = m->data;

	mData[0] *= s;
	mData[4] *= s;
	mData[8] *= s;
}


void pGslMatrixAddDiag(gsl_matrix* a, gsl_matrix* b)
{
	// a---> a+b.

	double* __restrict aData = a->data;
	double* __restrict bData = b->data;


	aData[0] += bData[0];
	aData[4] += bData[4];
	aData[8] += bData[8];
}


void adaptBackgroundThisModel(IplImage* ImAvoid, t_monster* pMons, double alpha, IplImage* I, IplImage* Ilhood, t_nrml aGaussian[])
{
	// Adapts mean for each pixel as:
	// mu = (1-alpha)*mu + I and similarly for beta, where beta comes from the bg likelihood image

	int V = I->height;
	int U = I->width;
	int Iws = I->widthStep;

	unsigned char* __restrict aRowI = (unsigned char*) I->imageData;
	float* __restrict aRowIlhood = (float*) Ilhood->imageData ;
	unsigned char* __restrict aRowImAvoid = (unsigned char*) ImAvoid->imageData;
	
	int wsImAvoid = ImAvoid->widthStep;
	int Vunrot = g_pSetts->V;
	double scFctr = ((double) Vunrot)/180.0;
	double minVariance = g_pSetts->minVariance * scFctr;


	double oneMinusAlpha = 1-alpha;

	int iPxl = 0;
	for(int v=0; v<V; ++v)
	{

		for(int u=0; u<U; ++u)
		{			

			if(aRowImAvoid[u] == (unsigned char) 0)
			{

				t_nrml* pGauss = &aGaussian[iPxl];

				// This pixel's value:
				int uTimes3 = u*3;
				double B = (double) aRowI[uTimes3++];
				double G = (double) aRowI[uTimes3++];
				double R = (double) aRowI[uTimes3];


				//------------------//
				// 1) Adapt mean	//
				//------------------//
			
				double* __restrict mData = pGauss->m->data;			
				double mB = mData[0];
				double mG = mData[1];
				double mR = mData[2];

				double adapted0 = (oneMinusAlpha)*mB + alpha*B;
				double adapted1 = (oneMinusAlpha)*mG + alpha*G;
				double adapted2 = (oneMinusAlpha)*mR + alpha*R;

				double beta = aRowIlhood[u] ;
				double oneMinusBeta = 1-beta;
				mData[0] = (oneMinusBeta)*mB + beta*adapted0;
				mData[1] = (oneMinusBeta)*mG + beta*adapted1;
				mData[2] = (oneMinusBeta)*mR + beta*adapted2;




				gsl_matrix* eMat = pMons->eMat;
				gsl_matrix* ee = pMons->ee;
				gsl_matrix* SLU = pMons->SLU;
				gsl_matrix* SLamda = pMons->SLamda;
				gsl_matrix* Stminus1 = pMons->Stminus1;
				gsl_permutation* p3 = pMons->p3;
				int s;
			
				gsl_matrix* S = pGauss->S;
				double& Sdet = pGauss->Sdet;
				gsl_matrix* Sinv = pGauss->Sinv;


				double* __restrict eMatData = eMat->data;

				eMatData[0] = B - mData[0] ;
				eMatData[1] = G - mData[1] ;
				eMatData[2] = R - mData[2] ;


			

#ifdef LF_FULL_COVAR_TOO
				gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, eMat, eMat, 0.0, ee);


				gsl_matrix_memcpy(Stminus1, S);


				gsl_matrix_memcpy(SLamda, Stminus1);
				gsl_matrix_scale(SLamda, oneMinusAlpha);
				gsl_matrix_scale(ee, alpha);
				gsl_matrix_add(SLamda, ee);



				gsl_matrix_scale(Stminus1, oneMinusBeta);
				gsl_matrix_scale(SLamda, beta);
				gsl_matrix_memcpy(S, Stminus1);
				gsl_matrix_add(S, SLamda);


			


			
//#ifdef LF_FULL_COVAR_TOO
				// TTT: minVariance = 16.0, config file
				S->data[0] = std::max(minVariance, S->data[0]);	
				S->data[4] = std::max(minVariance, S->data[4]);	
				S->data[8] = std::max(minVariance, S->data[8]);	


				gsl_matrix_memcpy(SLU, S);
				gsl_matrix_scale(SLU, 9); // TTT
				gsl_linalg_LU_decomp(SLU, p3, &s);
				Sdet = gsl_linalg_LU_det(SLU, s);
				gsl_linalg_LU_invert(SLU, p3, Sinv);
#endif



#ifdef LF_FULL_COVAR_TOO_1

				double* __restrict eeData = ee->data;
				eeData[0] = eMatData[0]*eMatData[0];
				eeData[4] = eMatData[1]*eMatData[1];
				eeData[8] = eMatData[2]*eMatData[2];


				pGslMatrixMemcpyDiag(Stminus1, S);


				pGslMatrixMemcpyDiag(SLamda, Stminus1);
				pGslMatrixScaleDiag(SLamda, oneMinusAlpha);
				pGslMatrixScaleDiag(ee, alpha);
				pGslMatrixAddDiag(SLamda, ee);



				pGslMatrixScaleDiag(Stminus1, oneMinusBeta);
				pGslMatrixScaleDiag(SLamda, beta);
				pGslMatrixMemcpyDiag(S, Stminus1);
				pGslMatrixAddDiag(S, SLamda);





				double* __restrict Sdata = S->data;

			
				// TTT: minVariance = 16.0, config file
				Sdata[0] = std::max(minVariance, Sdata[0]);	
				Sdata[4] = std::max(minVariance, Sdata[4]);	
				Sdata[8] = std::max(minVariance, Sdata[8]);	
			

			
				Sdet = 9*(Sdata[0] + Sdata[4] + Sdata[8]);

				double* __restrict SinvData = Sinv->data;
				SinvData[0] = 1/(9*Sdata[0]);
				SinvData[4] = 1/(9*Sdata[4]);
				SinvData[8] = 1/(9*Sdata[8]);		


#endif

				} //if ImAvoid==0
				
			++iPxl;
				
		
		} //u 

		aRowI += Iws;
		aRowIlhood += U;
		aRowImAvoid += wsImAvoid ;

	} //v



}



void showBackgroundAsImage(t_bg* pBg)
{
	int U = pBg->U;
	int V = pBg->V;
	
	IplImage* Im = pBg->Imean;
	int nChannels = Im->nChannels;
	float* aRowI = (float*) Im->imageData;
	IplImage* Isdet = pBg->Isdet;
	float* aRowIsdet = (float*) Isdet->imageData;
	
	int iPxl = 0;
	for(int v=0; v<V; ++v)
	{
		 

		for(int u=0; u<U; ++u)
		{
			double Sdet = pBg->aRGBdiag[iPxl].Sdet;

#ifdef LF_FULL_COVAR_TOO
			gsl_vector* m = pBg->aRGB[iPxl].m;
#endif
#ifdef LF_FULL_COVAR_TOO_1
			gsl_vector* m = pBg->aRGBdiag[iPxl].m;
#endif
			int uTimes3 = u*3;
			aRowI[uTimes3 + 0] = m->data[0]/255;
			aRowI[uTimes3 + 1] = m->data[1]/255;
			aRowI[uTimes3 + 2] = m->data[2]/255;

			aRowIsdet[u] = (float) Sdet;

			++iPxl;
		}

		aRowI += nChannels*U;
		aRowIsdet += U;
	}

	normalise32FImage(Isdet, 0, 1);

	//cvShowImage("2) Background Model", Im);
	//displayUnrotatedImage("2) Background Model", Im);
	//cvShowImage("Sdet", Isdet);
	//cvWaitKey(1);

}




