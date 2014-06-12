// Function declarations for poseidontest.
// Paris Kaimakis 19 Jan 2013


#ifndef _POSEIDONTEST_H_
#define _POSEIDONTEST_H_


#include "foreground.h"
#include "videocapture.h"
#include "datasetinfo.h"



int PoseidonTest(int firstArg, char* name);
void createVideoCapture(std::string szPath, std::string szFileName, const int U, const int V, t_vid* pVid);
void deleteVideoCapture(t_vid* pVid);
void loadVideo(t_vid* pVid);
void unloadVideo(t_vid* pVid);
int getNewFrame(t_vid* pVid);
int handleKeyInput(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg);
void processVideo(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track* pTrack);
void processVideoTestSet(std::string szPath, int U, int V, t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track* pTrack);
void createDataSetInfo(t_info* pInfo);
void displayInfo(t_info* pInfo);
void getPath(const std::string szPathWithExe, std::string& szPath);
int doEverything(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[]);
void initialiseEverythingParis(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg, t_track aTrack[]);



#endif //_POSEIDONTEST_H_


