// Monster struct implementation
// Paris Kaimakis 4 Jan 2013

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include "config.h"
#include "foreground.h"

extern t_settings* g_pSetts;

using std::cout;
using std::endl;

void createForegroundMonster(int U, int V, t_monster* pMons)
{
	pMons->eMat = gsl_matrix_calloc(3,1);
	pMons->ee = gsl_matrix_calloc(3,3);
	pMons->SLU = gsl_matrix_calloc(3,3);
	pMons->SLamda = gsl_matrix_calloc(3,3);
	pMons->Stminus1 = gsl_matrix_calloc(3,3);
	pMons->p3 = gsl_permutation_calloc(3);

	pMons->I32F1 = cvCreateImage(cvSize(U,V), IPL_DEPTH_32F, 1);
	pMons->I32F2 = cvCreateImage(cvSize(U,V), IPL_DEPTH_32F, 1);
	pMons->ImSea = cvCreateImage(cvSize(U,V), IPL_DEPTH_8U, 1);

	double scFct = ((double)V)/180;
	double sigma = 3.0*scFct; //3.0
	double gaussianMinValue = 0.1;
	int kerSizeHalf = (int) (sqrt( -2*sigma*sigma*std::log(gaussianMinValue) ) + 0.5);
	int kerSize = 2*kerSizeHalf + 1;

	pMons->sigma = sigma;
	pMons->kerSize = kerSize;
	pMons->IwithBorder = cvCreateImage(cvSize(U+kerSize-1,V+kerSize-1), IPL_DEPTH_32F, 1);


	int Vunrot = g_pSetts->V;
	double scFctr = ((double) Vunrot)/180.0;
	double sigmaWorstCase = g_pSetts->DeltaStDev/sqrt((double)Vunrot)*sqrt((double)V) *scFctr + g_pSetts->sigmaMin ;
	double thSigma = g_pSetts->thSigma;
	double kWorstCase = sqrt(-2*sigmaWorstCase*sigmaWorstCase*log(thSigma));
	

	int k = (int) (kWorstCase + 0.5);
	int K = 2*k+1;	
	pMons->k = k;				
	pMons->K = K;	// MUST be odd number.
	pMons->Ik = cvCreateImage(cvSize(K,K), IPL_DEPTH_32F, 1);
	pMons->IkTmp = cvCloneImage(pMons->Ik);
	
	int Uin = U+2*k;
	int Vin = V+2*k;
	pMons->IinForConvo = cvCreateImage(cvSize(Uin, Vin), IPL_DEPTH_32F, 1);

	pMons->aIk[0] = cvCreateImage(cvSize(K,1), IPL_DEPTH_32F, 1);
	pMons->aIk[1] = cvCreateImage(cvSize(1,K), IPL_DEPTH_32F, 1);
	pMons->aIin[0] = cvCreateImage(cvSize(Uin, V), IPL_DEPTH_32F, 1);
	pMons->aIin[1] = cvCreateImage(cvSize(U, Vin), IPL_DEPTH_32F, 1);


	pMons->I3chan32F = cvCreateImage(cvSize(U, V), IPL_DEPTH_32F, 3);
	pMons->I3chan8U = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 3);
	pMons->I3chan8U2 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 3);
	pMons->I1chan32F = cvCreateImage(cvSize(U, V), IPL_DEPTH_32F, 1);
	pMons->I1chan8U = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	pMons->I1chan8U2 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);
	pMons->I1chan8U3 = cvCreateImage(cvSize(U, V), IPL_DEPTH_8U, 1);

	int h=(int)((double)V/10 + .5);
	pMons->hx=h;
	pMons->hy=h;
	int hx = pMons->hx;
	int hy = pMons->hy;
	pMons->Itemplate = cvCreateImage(cvSize(hx,hy), IPL_DEPTH_32F, 3);
	pMons->ItempMatchResults = cvCreateImage(cvSize(U-hx+1,V-hy+1), IPL_DEPTH_32F, 1);
	
	pMons->Itemplate8U = cvCreateImage(cvSize(hx,hy), IPL_DEPTH_8U, 3);



	pMons->IorigClusterMask = cvCreateImage(cvSize(U+2,V+2), IPL_DEPTH_8U, 1);

}



void destroyForegroundMonster(t_monster* pMons)
{
	gsl_matrix_free(pMons->eMat);
	gsl_matrix_free(pMons->ee);
	gsl_matrix_free(pMons->SLU);
	gsl_matrix_free(pMons->SLamda);
	gsl_matrix_free(pMons->Stminus1);
	gsl_permutation_free(pMons->p3);

	cvReleaseImage(&pMons->I32F1);
	cvReleaseImage(&pMons->I32F2);
	cvReleaseImage(&pMons->ImSea);
	cvReleaseImage(&pMons->IwithBorder);


	cvReleaseImage(&pMons->Ik);
	cvReleaseImage(&pMons->IkTmp);
	cvReleaseImage(&pMons->IinForConvo);
	cvReleaseImage(&pMons->aIk[0]);
	cvReleaseImage(&pMons->aIk[1]);
	cvReleaseImage(&pMons->aIin[0]);
	cvReleaseImage(&pMons->aIin[1]);


	cvReleaseImage(&pMons->I3chan32F);
	cvReleaseImage(&pMons->I3chan8U);
	cvReleaseImage(&pMons->I3chan8U2);
	cvReleaseImage(&pMons->I1chan32F);
	cvReleaseImage(&pMons->I1chan8U);
	cvReleaseImage(&pMons->I1chan8U2);
	cvReleaseImage(&pMons->I1chan8U3);
	cvReleaseImage(&pMons->Itemplate);
	cvReleaseImage(&pMons->Itemplate8U);
	cvReleaseImage(&pMons->ItempMatchResults);
	cvReleaseImage(&pMons->IorigClusterMask);
	
}



