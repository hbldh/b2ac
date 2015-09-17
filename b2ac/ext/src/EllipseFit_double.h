#ifndef ELLIPSEFIT_DOUBLE_H
#define ELLIPSEFIT_DOUBLE_H

#include "EllipseFit.h"
#include "Eigenmethods_double.h"

EllipseDouble FitEllipse_double(int32_t *x, int32_t *y, int16_t numberOfPoints);

void calculateMatricesMAndT_double(int64_t *scatterMatrix, double *T, double *M);

void matrixInverseSymmetric3by3_double(int64_t *matrix, double *inverse);

void calculateTMatrix_double(int64_t *sm, double *s3i, double *T);

void calculateSecondTermOfM_double(int64_t *sm, double *T, double *M);

void calculateCompleteM_double(int64_t *sm, double *M);

void conicToGeneralForm_double(double *cc,
		EllipseDouble *results);

#endif // ELLIPSEFIT_DOUBLE_H