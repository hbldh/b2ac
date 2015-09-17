#ifndef ELLIPSEFIT_FLOAT_H
#define ELLIPSEFIT_FLOAT_H

#include "EllipseFit.h"
#include "Eigenmethods_float.h"

EllipseFloat FitEllipse_float(int32_t *x, int32_t *y, int16_t numberOfPoints);

void calculateMatricesMAndT_float(int64_t *scatterMatrix, float *T, float *M);

void matrixInverseSymmetric3by3_float(int64_t *matrix, float *inverse);

void calculateTMatrix_float(int64_t *sm, float *s3i, float *T);

void calculateSecondTermOfM_float(int64_t *sm, float *T, float *M);

void calculateCompleteM_float(int64_t *sm, float *M);

void conicToGeneralForm_float(float *cc,
		EllipseFloat *results);

#endif // ELLIPSEFIT_FLOAT_H
