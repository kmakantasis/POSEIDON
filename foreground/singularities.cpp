// Covariance matrix singularity analysis
// Paris Kaimakis 4 Jan 2013.


#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "foreground.h"


using std::cout;
using std::endl;




void removeSaturationInducedSingularities(t_monster* pMons, gsl_vector* m, double smallValue, double& Sdet, gsl_matrix* S)
{

	//gsl_vector* m = aGauss[iPxl].m;
	if(m->data[0]==0 || m->data[1]==0 || m->data[2]==0 || m->data[0]==255 || m->data[1]==255 || m->data[2]==255)
	{
		double aDiag[3] = {S->data[0], S->data[4], S->data[8]};

		//cout << "Sdet before = " << Sdet << endl;
		//cout << "Sbefore"<< endl;
		//printGslMatrix(S);

		// Find which one is zero:
		int iZero = -1;
		for(int ii=0; ii<3; ++ii)
			if(aDiag[ii]==0)
			{
				//iZero = ii;
				S->data[ii*4] = smallValue;
			}
		/*
		// Find the smallest non-zero:
		double minV = LARGE_VALUE;
		for(int ii=0; ii<3; ++ii)
		{
			if(ii!=iZero)
				if(aDiag[ii]<minV)
					minV = aDiag[ii];
		}

		// Replace the zero-diagonal entry of S with minV:
		S->data[iZero*4] = minV;
		//*/

		// Recalculate Sdet:
		gsl_matrix* SLU = pMons->SLU;
		gsl_permutation* p3 = pMons->p3;
		int s;
		gsl_matrix_memcpy(SLU, S);
		gsl_linalg_LU_decomp(SLU, p3, &s);
		Sdet = gsl_linalg_LU_det(SLU, s);
		
		//cout << "Sdet after = " << Sdet << endl;
		//cout << "Safter"<< endl;
		//printGslMatrix(S);

		//cvWaitKey();
	}


}
