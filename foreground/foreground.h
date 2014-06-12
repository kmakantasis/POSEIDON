// Definition of my background model class
//
// Paris Kaimakis 19/10/12



#ifndef _FOREGROUND_H_
#define _FOREGROUND_H_

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include "foregroundconstants.h"

#include <opencv2\video\background_segm.hpp>


struct s_horizon
{

	double xM, yM, thetaHdeg;
	CvMat* affMat;

	IplImage* ImLandInSea;

	int U, V;

	// Rotated image stuff:
	int Urot, Vrot;
	int xBoxLeft, yBoxTop;
	double xMrot, yMrot;
	CvMat* affMatRot;



};
typedef struct s_horizon t_horizon;




struct s_monster
{
	gsl_matrix* eMat ;
	gsl_matrix* ee;
	gsl_matrix* SLU ;
	gsl_permutation* p3 ;
	gsl_matrix* SLamda;
	gsl_matrix* Stminus1;


	IplImage* I32F1;
	IplImage* I32F2;
	IplImage* ImSea;

	// Stuff to do with improveBackgroundLikelihood():
	double sigma;
	int kerSize;
	IplImage* IwithBorder;  // Bigger image by kerSize-1 so as to allow for better convo.
							

	// Smart convo:
	int k;				// half-size for smart kernel.
	int K;				// kernel size (=2*k+1)	
	IplImage* Ik;		
	IplImage* IkTmp;
	IplImage* IinForConvo;

	IplImage* aIk[2];
	IplImage* aIin[2];


	// Horizon:
	IplImage* I3chan32F;
	IplImage* I3chan8U;
	IplImage* I3chan8U2;
	IplImage* I1chan32F;
	IplImage* I1chan8U;
	IplImage* I1chan8U2;
	IplImage* I1chan8U3;
	IplImage* Itemplate;
	IplImage* Itemplate8U;
	IplImage* ItempMatchResults;
	IplImage* IorigClusterMask ;

	int hx, hy;

};
typedef struct s_monster t_monster;





struct s_nrml 
{
	// Necessaries: m, S, Sinv, Sdet
	gsl_vector* m;
	gsl_matrix *S, *Sinv;
	double Sdet;

	// Usefuls for calculating likelihoods:
	gsl_vector* e;
	gsl_vector* SinvTimese;

};
typedef struct s_nrml t_nrml;



struct s_paper
{
	// For paper:
	
	double m;
	int t;

	IplImage* IpixelValuesOverTime;
	int Bbefore;
	int Gbefore;
	int Rbefore;

	// Visualise covariances:
	IplImage* IsamplesBG;
	IplImage* IsamplesGR;
	IplImage* IsamplesRB;


};
typedef struct s_paper t_paper;


struct s_track
{
	// The track taken from ONE target.

	double* aU;
	double* aV;

	IplImage* Itrack;	

};
typedef struct s_track t_track;



struct s_bg
{
	t_monster mons;
	IplImage* aIsample[NUM_FRAMES_FOR_TRAINING];
	int U, V;

	t_nrml* aRGB;
	t_nrml* aRGBdiag;

	IplImage* Imean;
	IplImage* Isdet;
	IplImage* Ilhood;
	IplImage* IlhoodDiag;
	IplImage* ImaskSdetZero;
	IplImage* ImaskRubbish;

	cv::BackgroundSubtractorMOG mog;
	cv::Mat mFrame, mForeground;


	IplImage* ImaskSea;
	IplImage* ImaskSky;
	IplImage* IforSegmo;
	IplImage* IsSignificant;
	IplImage* IlImproved2;


	IplConvKernel* kernel;


	// Important images:
	IplImage* Iresized;
	IplImage* Iresults;
	IplImage* ItoProcess;

	IplImage* ImLandInSea;

	t_paper paper;

	// Output to be used to form the posterior:
	IplImage* IjointBG;

};
typedef struct s_bg t_bg;








#endif // _FOREGROUND_H_



