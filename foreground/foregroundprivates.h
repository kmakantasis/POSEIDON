// Private functions available only to foreground project.
// Paris Kaimakis 17 Nov 2012

#ifndef _FOREGROUNDPRIVATES_H_
#define _FOREGROUNDPRIVATES_H_

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "foreground.h"


void recordBgTrainingSet(CvCapture* capture, IplImage* Iframe, t_horizon* pHorizon, IplImage* Iresized, CvMat* affMat, IplImage* ItoProcess, t_bg* pBg);
void calculateBgModel(int flag, t_bg* pBg);
void trainPixelMeans(int flag, t_bg* pBg);
void trainPixelCovars(int flag, t_bg* pBg);
void showBackgroundAsImage(t_bg* pBg);
void backgroundLikelihood(t_bg* pBg, IplImage* ImAvoid, IplImage* Iframe, IplImage* ImaskSky);
void adaptBackground(t_monster* pMons, double alpha, IplImage* I, t_bg* pBg);
void adaptBackgroundThisModel(IplImage* Isky, t_monster* pMons, double alpha, IplImage* I, IplImage* Ilhood, t_nrml aGaussian[]);
void makeHorizonHorizontal(IplImage* Isrc, CvMat* affMat, IplImage* Idst);
float pTotalKernelIntensity(int Keff, IplImage* I);
void pRecalculateKernel(int Keff, double sigma, IplImage* Ik);
float pixelColourLikelihoodFast(double B, double G, double R, t_nrml* gaussian);
void createNormal(t_nrml* pNrml, int dim);
void deleteNormal(t_nrml* pNrml);
void printGaussian(t_nrml* pNrml);
void createForegroundMonster(int U, int V, t_monster* pMons);
void destroyForegroundMonster(t_monster* pMons);
void removeSaturationInducedSingularities(t_monster* pMons, gsl_vector* m, double smallValue, double& Sdet, gsl_matrix* S);
int countNumberOfWhitePixels(IplImage* I);
void printGslMatrix(gsl_matrix* m);
void printGslVector(gsl_vector* v);
void normalise32FImage(IplImage* I, double minV, double maxV);
void maskToSilhouette(IplImage* I);
void equaliseImage(CvRect rect, IplImage* Ibgr, IplImage* IbgrEqual);
void linearTransform(IplImage* I, double alpha, double beta);
void pSmartSmoothNew(t_monster* pMons, IplImage* I, IplImage* Iout);
void pSmartSmoothNewFast(t_monster* pMons, double yMh, IplImage* ImAvoid, IplImage* I, IplImage* Iout);
void makeHorizonHorizontalPaddedWay(t_horizon* pHorizon, IplImage* I, IplImage* Ipadded, IplImage* Irot);
void displayUnrotatedImage(char* szText, IplImage* Irot);
void returnUnrotatedImage(char* szText, IplImage* Irot, IplImage* Iunrot);
void createStuffForPaper(t_paper* pPaper);
void deleteStuffForPaper(t_paper* pPaper);
void displayPixelIntensityGraphsForPaper(IplImage* Ibgr, int iFrame, t_paper* pPaper);
void tidyTrackArray(int nTargets, int iFrame, t_track aTrack[]);


#endif //_FOREGROUNDPRIVATES_H_

