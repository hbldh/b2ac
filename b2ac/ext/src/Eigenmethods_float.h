#ifndef EIGENMETHODS_FLOAT_H_INCLUDED
#define EIGENMETHODS_FLOAT_H_INCLUDED

#include <math.h>
#include <stdint.h>

#define QR_ALGORITHM_TOLERANCE_FLOAT (1e-4)

void QRAlgorithm_float(float *A, float *eigenvalues);

void matrixToHessenbergForm_Givens_float(float *A);

void matrixToHessenbergForm_Householder_float(float *A);

void givensRotation_float(float *csVect, float a, float b);

void performQRStep_Givens_float(float *A, int m);

void inverseIteration_float(float *A, float* eigenvector,
		float eigenvalue, int32_t numberOfIterations);

float matrixInverse3by3_float(float *A, float *inverse);

#endif // EIGENMETHODS_FLOAT_H_INCLUDED
