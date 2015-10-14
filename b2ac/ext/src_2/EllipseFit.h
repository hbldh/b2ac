#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

typedef struct _EllipseDouble
{
    double centerX;
    double centerY;
    double majorAxis;
    double minorAxis;
    double rotationAngle;
	bool isValid;
} EllipseDouble;

EllipseDouble FitEllipse_double(int32_t *x, int32_t *y, int16_t numberOfPoints);

EllipseDouble FitEllipse_int(int32_t *x, int32_t *y, int16_t numberOfPoints);

int32_t removeMeanFromCoordinates(int32_t *coords, int16_t numberOfPoints);

void calculateScatterMatrix(int32_t *x, int32_t *y, int16_t numberOfPoints, int32_t *scatterMatrix);

void matrixInverseSymmetric3by3_double(int32_t *matrix, double *inverse);

int64_t matrixInverseSymmetric3by3_int(int32_t *matrix, int64_t *inverse);

void calculateTMatrix_double(int32_t *sm, double *s3i, double *T);

void calculateTMatrix_int(int32_t *sm, int64_t *s3i, int64_t *T);

void calculateSecondTermOfM_double(int32_t *sm, double *T, double *M);

void calculateSecondTermOfM_int(int32_t *sm, int64_t *T, int32_t *M, int64_t determinant);

void calculateCompleteM_double(int32_t *sm, double *M);

void calculateCompleteM_int(int32_t *sm, int32_t *M);

void conicToGeneralForm_double(double *conicCoefficients, EllipseDouble *results);


