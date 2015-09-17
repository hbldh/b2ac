#ifndef ELLIPSEFIT_INT_H
#define ELLIPSEFIT_INT_H

#include <stdlib.h>
#include "EllipseFit.h"
#include "Eigenmethods_int.h"

EllipseFloat FitEllipse_int(int32_t *x,
							   int32_t *y,
							   int16_t numberOfPoints);

int64_t calculateMatricesMAndT_int(int64_t *scatterMatrix,
									  int64_t *T,
									  int64_t *M);

int64_t matrixInverseSymmetric3by3_int(int64_t *A,
										   int64_t *inverse);

void calculateTMatrix_int(int64_t *scatterMatrix,
							 int64_t *s3Inverse,
							 int64_t *T);

int64_t scaleTMatrix(int64_t *T, int64_t S3Determinant);

void calculateSecondTermOfM_int(int64_t *scatterMatrix,
									int64_t *T,
									int64_t *M,
									int64_t determinant);

void calculateCompleteM_int(int64_t *scatterMatrix, int64_t *M);

void conicToGeneralForm_int(int64_t *cc, EllipseFloat *results);

float convertSqrtToFloat(uint64_t val);

#endif // ELLIPSEFIT_INT_H
