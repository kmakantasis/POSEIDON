// Poseidon main function
//
// Paris Kaimakis 12 Oct 2012

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>
#include "config.h"
#include "foreground.h"
#include "foregroundinterface.h"
#include "videocapture.h"
#include "poseidontest.h"


using std::cout;
using std::endl;


#ifdef NO_INTEGRATION
t_settings* g_pSetts;
t_horizon* g_pHorizon;


int main1(int argc, char* argv[])
{

	char* name = argv[1];
	int firstArg = argc;

	t_info info, *pInfo=&info;
	t_vid vid, *pVid=&vid;
	t_bg bg, *pBg=&bg;
	t_horizon horizon, *pHorizon=&horizon;
	g_pHorizon = pHorizon;
	t_track aTrack[MAX_TARGETS];
	t_settings settings;
	g_pSetts=&settings;
	std::string szPath;

	getPath(name, szPath);
	readConfigFile(szPath, g_pSetts);
	createDataSetInfo(pInfo);
	createHorizon(g_pSetts->U, g_pSetts->V, pHorizon);


	if(firstArg==1)
		processVideoTestSet(szPath, g_pSetts->U, g_pSetts->V, pInfo, pVid, pHorizon, pBg, aTrack);
	if(firstArg==2)
	{
		std::string szFilename = name;
		std::string szFilenameWithPath = szPath + szFilename ;
		createVideoCapture(szPath, szFilenameWithPath, g_pSetts->U, g_pSetts->V, pVid);

		// Create Bg class, Process Video:
		processVideo(pInfo, pVid, pHorizon, pBg, aTrack);


		deleteVideoCapture(pVid);
	}


	deleteHorizon(pHorizon);


	
//*/
	return 0;

}

#endif

