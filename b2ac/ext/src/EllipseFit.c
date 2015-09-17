#include "EllipseFit.h"

int32_t findMeanOfCoordinates(int32_t *coordinates, int16_t numberOfPoints)
{
    int16_t _i;
    int32_t _meanPoint = 0;

    for(_i = 0; _i < numberOfPoints; _i++)
    {
        _meanPoint += coordinates[_i];
    }

    return (_meanPoint / numberOfPoints);
}

void calculateScatterMatrix(int32_t *x, int32_t *y, int32_t xMean,
		int32_t yMean, int16_t numberOfPoints, int64_t *scatterMatrix)
{
    int16_t _i;

    int64_t _x;
    int64_t _y;

    int64_t tmpX2;
    int64_t tmpX3;
    int64_t tmpY2;
    int64_t tmpY3;

    for(_i = 0; _i < numberOfPoints; _i++)
    {
    	_x = (int64_t) (x[_i] - xMean);
    	_y = (int64_t) (y[_i] - yMean);

        tmpX2 = _x * _x;
        tmpX3 = tmpX2 * _x;
        tmpY2 = _y * _y;
        tmpY3 = tmpY2 * _y;

        // First row of scatter matrix
        scatterMatrix[0] += tmpX2 * tmpX2;
        scatterMatrix[1] += tmpX3 * _y;
        scatterMatrix[2] += tmpX2 * tmpY2;
        scatterMatrix[3] += tmpX3;
        scatterMatrix[4] += tmpX2 * _y;
        scatterMatrix[5] += tmpX2;

        // Second row of scatter matrix
        // Index 6 is a duplicate of index 2.
        scatterMatrix[7] += tmpY3 * _x;
        // Index 8 is a duplicate of index 4.
        scatterMatrix[9] += tmpY2 * _x;
        scatterMatrix[10] += _x * _y;

        // Third row of scatter matrix
        scatterMatrix[11] += tmpY2 * tmpY2;
        // Index 12 is a duplicate of index 9.
        scatterMatrix[13] += tmpY2 * _y;
        scatterMatrix[14] += tmpY2;

        // Fourth row of scatter matrix
        // Index 15 is a duplicate of index 5.
        // Index 16 is a duplicate of index 10.
        scatterMatrix[17] += _x;

        // Fifth row of scatter matrix
        // Index 18 is a duplicate of index 14.
        scatterMatrix[19] += _y;
    }

    // Set duplicates.
    scatterMatrix[6] = scatterMatrix[2];
    scatterMatrix[8] = scatterMatrix[4];
    scatterMatrix[12] = scatterMatrix[9];
    scatterMatrix[15] = scatterMatrix[5];
    scatterMatrix[16] = scatterMatrix[10];
    scatterMatrix[18] = scatterMatrix[14];

    // Set number of points in last slot.
    scatterMatrix[20] = (int64_t) numberOfPoints;
}

