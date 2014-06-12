// Stuff to do with Gaussian distribo for each pixel
// Paris Kaimakis 26 Oct 2012

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "foreground.h"
#include "foregroundprivates.h"



using std::cout;
using std::endl;


void createNormal(t_nrml* pNrml, int dim)
{

	pNrml->m = gsl_vector_calloc(dim);
	pNrml->S = gsl_matrix_calloc(dim, dim);
	pNrml->Sinv = gsl_matrix_calloc(dim, dim);
	pNrml->Sdet = 0;

	pNrml->e = gsl_vector_calloc(dim);
	pNrml->SinvTimese = gsl_vector_calloc(dim);


}



void deleteNormal(t_nrml* pNrml)
{
	gsl_vector_free(pNrml->m);
	gsl_matrix_free(pNrml->S);
	gsl_matrix_free(pNrml->Sinv);

	gsl_vector_free(pNrml->e);
	gsl_vector_free(pNrml->SinvTimese);
}



void printGaussian(t_nrml* pNrml)
{


	cout << "m = " << endl;
	printGslVector(pNrml->m);

	cout << "S = " << endl;
	printGslMatrix(pNrml->S);

	cout << "Sinv = " << endl;
	printGslMatrix(pNrml->Sinv);

	cout << "Sdet = " << pNrml->Sdet << endl << endl;


}






