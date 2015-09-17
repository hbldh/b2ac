#ifndef ELLIPSEFIT_H
#define ELLIPSEFIT_H

#include <stdint.h>
#include <math.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
// Do not add anything in Windows.
#else
#include <stdbool.h>
#endif


typedef struct _EllipseDouble
{
    double centerX;
    double centerY;
    double xAxis;
    double yAxis;
    double rotationAngle;
	bool isValid;
} EllipseDouble;

typedef struct _EllipseFloat
{
    float centerX;
    float centerY;
    float xAxis;
    float yAxis;
    float rotationAngle;
	bool isValid;
} EllipseFloat;

int32_t findMeanOfCoordinates(int32_t *coordinates, int16_t numberOfPoints);

void calculateScatterMatrix(int32_t *x, int32_t *y, int32_t xMean,
		int32_t yMean, int16_t numberOfPoints, int64_t *scatterMatrix);

#endif // ELLIPSEFIT_H

