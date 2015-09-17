#ifndef EIGENMETHODS_DOUBLE_H_INCLUDED
#define EIGENMETHODS_DOUBLE_H_INCLUDED

#include <math.h>
#include <stdint.h>

#define QR_ALGORITHM_TOLERANCE_DOUBLE (1e-10)

void QRAlgorithm_double(double *A, double *eigenvalues);

void matrixToHessenbergForm_Givens_double(double *A);

void matrixToHessenbergForm_Householder_double(double *A);

void givensRotation_double(double *csVect, double a, double b);

void performQRStep_Givens_double(double *A, int m);

void inverseIteration_double(double *A, double* eigenvector,
		double eigenvalue, int32_t numberOfIterations);

double matrixInverse3by3_double(double *A, double *inverse);

#endif // EIGENMETHODS_DOUBLE_H_INCLUDED
