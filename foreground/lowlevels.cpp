// Low level stuff here.
// Paris Kaimakis 7 Jan 2013


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


using std::cout;
using std::endl;


int countNumberOfWhitePixels(IplImage* I)
{
	// I is single-channel, 8-bit.


	int nPixels=0;
	int V = I->height;
	int U = I->width;
	unsigned char* aRow = (unsigned char*) I->imageData;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRow[u]== (unsigned char) 255)
				++nPixels;
		} //u

		aRow += U;

	} //v

	return nPixels;
}


void printGslVector(gsl_vector* v)
{
	for(int i=0; i<v->size; ++i)
		cout << "v[" << i << "] = " << gsl_vector_get(v, i) << "\t" << endl;

	cout << endl << endl;
}



void printGslMatrix(gsl_matrix* m)
{
	for(int i=0; i<m->size1; ++i)
	{
		for(int j=0; j<m->size2; ++j)
			cout << "m[" << i << "," << j << "] = " << gsl_matrix_get(m, i, j) << "\t";

		cout << endl;
	}

	cout << endl << endl;
}




void normalise32FImage(IplImage* I, double minV, double maxV)
{
	double minVnow, maxVnow;
	cvMinMaxLoc(I, &minVnow, &maxVnow);
	cvScale(I, I, 1.0, -minVnow);
	cvScale(I, I, (maxV-minV)/(maxVnow-minVnow), 0.0);
	cvScale(I, I, 1.0, minV);


}



void maskToSilhouette(IplImage* I)
{

	int V = I->height;
	int U = I->width;
	int ws = I->widthStep;

	unsigned char* aRow = (unsigned char*) I->imageData;

	for(int v=0; v<V; ++v)
	{
		for(int u=0; u<U; ++u)
		{
			if(aRow[u]== (unsigned char) 255)
				aRow[u] = (unsigned char) 0;
			else
				aRow[u] = (unsigned char) 255;

		}

		aRow += ws;
	}
}


