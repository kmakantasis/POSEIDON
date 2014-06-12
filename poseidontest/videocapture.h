// Everything to do with videocapture for poseidontest
// Paris Kaimakis 15/11/12


#ifndef _VIDEOCAPTURE_H_
#define _VIDEOCAPTURE_H_

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>




#define NUM_TESTSET_VIDEOS	10


#define M_FILE				0	// obsolete
#define M_TESTSET			1	// obsolete
#define M_WEB				2


struct s_vid
{
	std::string szFileNameData;		// Filename with path, or "", or "web".
	std::string szPath;				

	CvCapture* capture;
	IplImage* Iframe;

	int iFrame;

	IplImage* Iresized;
	IplImage* Iresults;
	IplImage* ItoProcess;

	CvVideoWriter *writer;

	int tMode; // eg wem-mode

};
typedef struct s_vid t_vid;




#endif //_VIDEOCAPTURE_H_

