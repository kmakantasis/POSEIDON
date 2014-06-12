// Key input
// Paris Kaimakis 26 Oct 2012


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "foreground.h"
#include "foregroundinterface.h"
#include "videocapture.h"
#include "datasetinfo.h"
#include "poseidontest.h"

using std::cout;
using std::endl;




int handleKeyInput(t_info* pInfo, t_vid* pVid, t_horizon* pHorizon, t_bg* pBg)
{
	int c = cvWaitKey(1);
	int bContinue = 1;

	if(c==-1)
		return bContinue;


	switch(c)
	{
	case 27:
		// Escape:
		bContinue = 0;
		break;
	case 'i':
	case 'I':
		displayInfo(pInfo);
		cout << endl << endl;
		cout << "Please select option" << endl;
		break;
	case 'q':
	case 'Q':
		cvWaitKey();
		break;
	case 't':
	case 'T':
		trainBgModel(pVid->capture, pVid->Iframe, pHorizon, pVid->Iresized, pHorizon->affMat, pVid->ItoProcess, pBg);
		break;
	case 'r':
	case 'R':
		pVid->iFrame = 0;
		cvSetCaptureProperty(pVid->capture, CV_CAP_PROP_POS_FRAMES, 1);
		break;
	}

	return bContinue;
}




