// Stuff to do with videocapture class
// Paris Kaimakis 26 Oct 2012


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "foregroundinterface.h"
#include "foregroundconstants.h"
#include "videocapture.h"
#include "poseidontest.h"


using std::cout;
using std::endl;

#define CAPTURE_FROM_AVI
#define F_RECORD_RESULTS




void createVideoCapture(std::string szPath, std::string szFileNameWithPath, const int U, const int V, t_vid* pVid)
{
	cout << "\tCreating video class... " ;
	cout.flush();
	
	pVid->szPath = szPath;
	


#ifdef CAPTURE_FROM_AVI
	pVid->szFileNameData = szFileNameWithPath;

	// (Re)extract the input argument:
	std::string szArg = "";
	int pos = szFileNameWithPath.rfind('\\');
	int length = szFileNameWithPath.length();
	for(int i=pos+1; i<length; ++i)
	{
		szArg += szFileNameWithPath[i];
	}

	
	if(szArg=="web")
	{
		pVid->tMode = M_WEB;

		pVid->szFileNameData = "http://root:poseidon@147.27.11.212:80/mjpg/video.mjpg?resolution=320x240&req_fps=25&.mjpg";
				// works, but not always

		//pVid->szFileNameData = "http://root:poseidon@147.27.11.212:80//axis-cgi/mjpg/video.cgi?resolution=640x480&req_fps=30&.mjpg";
		//pVid->szFileNameData = "http://root:poseidon@147.27.11.212:80/mjpg/video.mjpg?resolution=320x240&req_fps=20&.mjpg";
		//pVid->szFileNameData = "http://root:poseidon@147.27.11.212:80/axis-cgi/mjpg/video.cgi?resolution=320x240";

		// Tests -- All failed:
		//1) 
		//pVid->szFileNameData = "axrtpu://root:poseidon@147.27.11.212:80/mpeg4/media.amp";
		//2) 
		//pVid->szFileNameData = "axrtsp://root:poseidon@147.27.11.212:80/mpeg4/media.amp";
		// 3)
		//pVid->szFileNameData = "axrtsphttp://root:poseidon@147.27.11.212:80/axis-media/media.amp";
		// 4)
		//pVid->szFileNameData = "axrtpm://root:poseidon@147.27.11.212:80/mpeg4/media.amp";
		// 5) 
		//pVid->szFileNameData = "rtsp://root:poseidon@147.27.11.212:80/mpeg4/media.amp";
	}								
#else
	pVid->szFileNameData = "camera";
#endif


	CvSize size = cvSize(U, V);
	pVid->Iresized = cvCreateImage(size, IPL_DEPTH_8U, 3);
	pVid->Iresults = cvCreateImage(size, IPL_DEPTH_8U, 3);
	pVid->ItoProcess = cvCreateImage(size, IPL_DEPTH_8U, 3);

	
	cout << "Done! " << endl;

}




void deleteVideoCapture(t_vid* pVid)
{
	cout << "\tDestroying video class... " ;
	cout.flush();

	cvReleaseImage(&pVid->Iresized);
	cvReleaseImage(&pVid->Iresults);
	cvReleaseImage(&pVid->ItoProcess);

#ifdef F_RECORD_RESULTS
	cvReleaseVideoWriter(&pVid->writer);
#endif


	cout << "Done!" << endl;
}




void loadVideo(t_vid* pVid)
{
	// Loads the video indicated by pVid->szFileNameData.

	cout << "\tLoading video... " ;
	//cout.flush();

#ifdef CAPTURE_FROM_AVI
	pVid->capture = cvCaptureFromAVI(pVid->szFileNameData.c_str());
#else
	pVid->capture = cvCaptureFromCAM(0);
#endif


	// Recorder:
#ifdef F_RECORD_RESULTS
	std::string szFileNameResults = pVid->szFileNameData;

	int n=szFileNameResults.length();
	szFileNameResults.insert(n-4, "results");
	double fps = cvGetCaptureProperty(pVid->capture, CV_CAP_PROP_FPS);	
	
	if(pVid->tMode==M_WEB) szFileNameResults = pVid->szPath + "webresults.avi";


	int U = pVid->Iresized->width;
	int V = pVid->Iresized->height;
	CvSize size = cvSize(U, V);

	//----------------------------------------------------------------------//
	// Format Wars!!														//
	// OCV = ability to open recorded video through OpenCV in retrospect	//
	// MPC = ability to play with Media Player Classic						//
	//																		//
	// Codec names from http://www.cs.iit.edu/~agam/cs512/lect-notes/opencv-intro/opencv-intro.html
	//----------------------------------------------------------------------//

	// MPEG-1: OCV+, MPC-, VLC+
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('P','I','M','1'), fps, size); // low quality
	
	// Motion-JPEG: OCV+ MPC+ Perhaps not great compression??
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('M','J','P','G'), fps, size);

	// MPEG-4: OCV-, MPC+
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('D','I','V','X'), fps, size);

	// MPEG-4.2: OCV+, MPC+
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('M','P','4','2'), fps, size);
					
	// MPEG-4.3: OCV+, MPC+
	pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('D','I','V','3'), fps, size);

	// H263: OCV-, MPC+
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('U','2','6','3'), fps, size);
	
	// H263I: OCV-, MPC-, VLC-, Output file is empty.
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('I','2','6','3'), fps, size);
	
	// FLV1: OCV+, MPC+, Perhaps best compression ratio
	// I think it gave me trouble ...?
	//pVid->writer = cvCreateVideoWriter(szFileNameResults.c_str(), CV_FOURCC('F','L','V','1'), fps, size);



#endif


	pVid->Iframe = cvQueryFrame(pVid->capture);
	pVid->iFrame = 0 + NUM_FRAMES_FOR_TRAINING;

	int Uraw = pVid->Iframe->width;
	int Vraw = pVid->Iframe->height;

	cout << "Done!" << endl;

}



void unloadVideo(t_vid* pVid)
{
	cvReleaseCapture(&pVid->capture);
}





int getNewFrame(t_vid* pVid)
{
	IplImage* Iframe = pVid->Iframe;
	Iframe = cvQueryFrame(pVid->capture);

	if(!Iframe)
		return 0;

	// Prosoxi!! To pio katw kanonika xreiazetai, alla to kanoume bypass logw tou integration
	cvResize(Iframe, pVid->Iresized);
	
	
	cvCopy(pVid->Iresized, pVid->Iresults);


	
	return 1;
}




