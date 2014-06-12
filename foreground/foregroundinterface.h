// Interface for the foreground object.
// Paris Kaimakis 17 Nov 2012

#ifndef _FOREGROUNDINTERFACE_H_
#define _FOREGROUNDINTERFACE_H_

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "config.h"
#include "foreground.h"



void createBg(t_bg* pBg, int U, int V);
void destroyBg(t_bg* pBg);
void trainBgModel(CvCapture* capture, IplImage* Iframe, t_horizon* pHorizon, IplImage* Iresized, CvMat* affMat, IplImage* ItoProcess, t_bg* pBg);
void trainBgModelNew(t_bg* pBg);
void showBackgroundAsImage(t_bg* pBg);
void adaptBackground(t_monster* pMons, double alpha, IplImage* I, t_bg* pBg);
void adaptBackgroundThisModel(IplImage* Isky, t_monster* pMons, double alpha, IplImage* I, IplImage* Ilhood, t_nrml aGaussian[]);
void removeBackground(t_monster* pMons, IplImage* ImLandInSea, double yMh, t_nrml aGaussian[], IplImage* Iframe, IplImage* Ilhood, IplImage* ImaskSky, IplImage* IforSegmo, IplImage* IjointBG);
void makeHorizonHorizontal(IplImage* Isrc, CvMat* affMat, IplImage* Idst);
void displayPixelIntensityGraphsForPaper(IplImage* Ibgr, int iFrame, t_paper* pPaper);
int postprocess(t_horizon* pHorizon, IplImage* Ilhood, IplImage* ImLandInSea, IplImage* ItoProcess, t_monster* pMons, IplConvKernel* kernel, IplImage* IlBack32F, IplImage* Iresized, IplImage* Isignif, int iFrame, t_track aTrack[]);
void postprocessTest(t_bg* pBg, IplImage* Ibgr);
void equaliseImage(CvRect rect, IplImage* Ibgr, IplImage* IbgrEqual);
void calculateImageStatistics(IplImage* I);
void linearTransform(IplImage* I, double alpha, double beta);
void finalAdaptations(IplImage* IsSignif, IplImage* Ilhood, IplImage* IlImproved2);
void readConfigFile(const std::string szPath, t_settings* pSets);
void createHorizon(int U, int V, t_horizon* pHorizon);
void deleteHorizon(t_horizon* pHorizon);
void findHorizon(CvCapture* pCapt, IplImage* Iraw, IplImage* I, t_horizon* pHorizon);
void makeHorizonHorizontalPaddedWay(t_horizon* pHorizon, IplImage* I, IplImage* Ipadded, IplImage* Irot);
void displayUnrotatedImage(char* szText, IplImage* Irot);
void returnUnrotatedImage(char* szText, IplImage* Irot, IplImage* Iunrot);
void createTrackArray(CvCapture* pCapture, int U, int V, t_track aTrack[]);
void deleteTrackArray(t_track aTrack[]);
void tidyTrackArray(int nTargets, int iFrame, t_track aTrack[]);
void drawTrackArrayFrameByFrame(int iFrame, t_track aTrack[]);
void drawTrackArrayOneOff(IplImage* Ibgr, int nFrames, t_track aTrack[]);


#endif //_FOREGROUNDINTERFACE_H_


