#include "stdafx.h"
#include "EllipseFit.h"

EllipseDouble FitEllipse_double(int32_t *x, int32_t *y, int16_t numberOfPoints)
{
    static int16_t _i;
    // Storage for the scatter matrix. It is a 6x6 symmetrical matrix.
    static int32_t *scatterMatrix = (int32_t *) malloc(sizeof(int32_t) * 21);
	// Storage for matrices to calculate on the way to eigensystem matrix.
    static double *S3Inverse = (double *) malloc(sizeof(double) * 6);
    static double *T = (double *) malloc(sizeof(double) * 9);
	// Eigensystem matrix storage.
	static double *M = (double *) malloc(sizeof(double) * 9);
	// Eigenvectors and eigenvalues storage.
	//static double *eigenvalues = (double *) malloc(sizeof(double) * 3);
	//static double *eigenvectors = (double *) malloc(sizeof(double) * 9);
	static double eigenvalues[3] = {1432841.99428,5476.75547089,-1122751.7107};
	static double eigenvectors[9] = {0.889048686155,-0.457588868649,0.00503636728127,-0.0143632566829,-0.00589873101016,-0.999969859516,-0.457587292769,-0.889144325889,-0.00590889697334};
	static double *conicCoefficients = (double *) malloc(sizeof(double) * 6);

    static int32_t _xMean = 0;
    static int32_t _yMean = 0;
	static int16_t ellipticalEigenvectorIndex = -1;

	static EllipseDouble results;

    // Reset scatter matrix. All other matrices will be replaced "in transit".
    for(_i = 0; _i < 21; _i++)
    {
        scatterMatrix[_i] = 0;
    }

    _xMean = removeMeanFromCoordinates(x, numberOfPoints);
    _yMean = removeMeanFromCoordinates(y, numberOfPoints);

	// Calculate the T and M matrices needed.
    calculateScatterMatrix(x, y, numberOfPoints, scatterMatrix);
    matrixInverseSymmetric3by3_double(&scatterMatrix[15], S3Inverse);
    calculateTMatrix_double(scatterMatrix, S3Inverse, T);
	calculateSecondTermOfM_double(scatterMatrix, T, M);
	calculateCompleteM_double(scatterMatrix, M);

	// Now we have M, apply eigensolver to get a1 vector.
	printf("Double version:\n");
	for(_i = 0; _i < 9; _i++)
    {
        printf("M[%d] = %f\n", _i, M[_i]);
    }
	printf("\n");

	// Eigenvalues calculated; evaluate which one is elliptical.
	// It can be proven that only one eigenvalue will be positive in case
	// of elliptical solution; if none fulfills this, then fitting was unsucessful.
	ellipticalEigenvectorIndex = -1;
	for(_i = 0; _i < 3; _i++)
	{
		if ((4 * eigenvectors[_i] * eigenvectors[6 + _i]) - (eigenvectors[3 + _i] * eigenvectors[3 + _i]) > 0.0)
		{
			ellipticalEigenvectorIndex = _i;
			break;
		}
	}

	if (ellipticalEigenvectorIndex == -1)
	{
		// Non-elliptical solution found. Returning empty results.
		results.isValid = false;
	}
	else
	{
		// Sort the coefficients calculated and find the three last ones.
		conicCoefficients[0] = eigenvectors[_i];
		conicCoefficients[1] = eigenvectors[3 +_i];
		conicCoefficients[2] = eigenvectors[6 + _i];
		conicCoefficients[3] = conicCoefficients[0] * T[0] + conicCoefficients[1] * T[1] + conicCoefficients[2] * T[2];
		conicCoefficients[4] = conicCoefficients[0] * T[3] + conicCoefficients[1] * T[4] + conicCoefficients[2] * T[5];
		conicCoefficients[5] = conicCoefficients[0] * T[6] + conicCoefficients[1] * T[7] + conicCoefficients[2] * T[8];

		printf("Double version:\n");
		for(_i = 0; _i < 6; _i++)
		{
			printf("c[%d] = %f\n", _i, conicCoefficients[_i]);
		}
		printf("\n");

		// Convert to general form and store this information
		// in the output struct.
		conicToGeneralForm_double(conicCoefficients, &results);
		// Add the previously subtracted mean point.
		results.centerX += _xMean;
		results.centerY += _yMean;
	}

	printf("Double version:\n");
	printf("Center = %f, %f\n", results.centerX, results.centerY);
	printf("Axes (Maj, Min) = %f, %f\n", results.majorAxis, results.minorAxis);
	printf("Rotation Angle = %f\n", results.rotationAngle);
	printf("Valid = %c\n", results.isValid);

	return results;
}

EllipseDouble FitEllipse_int(int32_t *x, int32_t *y, int16_t numberOfPoints)
{
    static int16_t _i;
    // Storage for the scatter matrix. It is a 6x6 symmetrical matrix.
    static int32_t *scatterMatrix = (int32_t *) malloc(sizeof(int32_t) * 21);
	// Storage for variables and matrices to calculate on the way to eigensystem matrix.
    static int64_t *S3Inverse = (int64_t *) malloc(sizeof(int64_t) * 6);
    static int64_t S3Determinant = 0;
	static int64_t *T = (int64_t *) malloc(sizeof(int64_t) * 9);
	// Eigensystem matrix storage.
	static int32_t *M = (int32_t *) malloc(sizeof(int32_t) * 9);
	// Eigenvectors and eigenvalues storage.

    static int32_t _xMean = 0;
    static int32_t _yMean = 0;

	static EllipseDouble results;

    // Reset scatter matrix. All other matrices will be replaced "in transit".
    for(_i = 0; _i < 21; _i++)
    {
        scatterMatrix[_i] = 0;
    }

    _xMean = removeMeanFromCoordinates(x, numberOfPoints);
    _yMean = removeMeanFromCoordinates(y, numberOfPoints);

	// Calculate the T and M matrices needed.
    calculateScatterMatrix(x, y, numberOfPoints, scatterMatrix);
    S3Determinant = matrixInverseSymmetric3by3_int(&scatterMatrix[15], S3Inverse);
    calculateTMatrix_int(scatterMatrix, S3Inverse, T);
	calculateSecondTermOfM_int(scatterMatrix, T, M, S3Determinant);
	calculateCompleteM_int(scatterMatrix, M);

	// Now we have M, apply eigensolver to get a1 vector.
	printf("Integer version:\n");
	for(_i = 0; _i < 9; _i++)
    {
        printf("M[%d] = %d\n", _i, M[_i]);
    }
	printf("\n");

	return results;
}

int32_t removeMeanFromCoordinates(int32_t *coords, int16_t numberOfPoints)
{
	return 0;
    static int16_t _i;
    static int32_t _meanPoint;

    _meanPoint = 0;
    for(_i = 0; _i < numberOfPoints; _i++)
    {
        _meanPoint += coords[_i];
    }

    _meanPoint /= numberOfPoints;

    for(_i = 0; _i < numberOfPoints; _i++)
    {
        coords[_i] -= _meanPoint;
    }

    return (_meanPoint);
}

void calculateScatterMatrix(int32_t *x, int32_t *y, int16_t numberOfPoints, int32_t *scatterMatrix)
{
    static int16_t _i;
    static int32_t tmpX2;
    static int32_t tmpX3;
    static int32_t tmpY2;
    static int32_t tmpY3;

    for(_i = 0; _i < numberOfPoints; _i++)
    {
        tmpX2 = x[_i] * x[_i];
        tmpX3 = tmpX2 * x[_i];
        tmpY2 = y[_i] * y[_i];
        tmpY3 = tmpY2 * y[_i];

        // First row of scatter matrix
        scatterMatrix[0] += tmpX2 * tmpX2;
        scatterMatrix[1] += tmpX3 * y[_i];
        scatterMatrix[2] += tmpX2 * tmpY2;
        scatterMatrix[3] += tmpX3;
        scatterMatrix[4] += tmpX2 * y[_i];
        scatterMatrix[5] += tmpX2;

        // Second row of scatter matrix
        // Index 6 is a duplicate of index 2.
        scatterMatrix[7] += tmpY3 * x[_i];
        // Index 8 is a duplicate of index 4.
        scatterMatrix[9] += tmpY2 * x[_i];
        scatterMatrix[10] += x[_i] * y[_i];

        // Third row of scatter matrix
        scatterMatrix[11] += tmpY2 * tmpY2;
        // Index 12 is a duplicate of index 9.
        scatterMatrix[13] += tmpY2 * y[_i];
        scatterMatrix[14] += tmpY2;

        // Fourth row of scatter matrix
        // Index 15 is a duplicate of index 5.
        // Index 16 is a duplicate of index 10.
        scatterMatrix[17] += x[_i];

        // Fifth row of scatter matrix
        // Index 18 is a duplicate of index 14.
        scatterMatrix[19] += y[_i];
    }

    // Set duplicates.
    scatterMatrix[6] = scatterMatrix[2];
    scatterMatrix[8] = scatterMatrix[4];
    scatterMatrix[12] = scatterMatrix[9];
    scatterMatrix[15] = scatterMatrix[5];
    scatterMatrix[16] = scatterMatrix[10];
    scatterMatrix[18] = scatterMatrix[14];

    // Set number of points in last slot.
    scatterMatrix[20] = (int32_t) numberOfPoints;
}

void matrixInverseSymmetric3by3_double(int32_t *matrix, double *inverse)
{
    static double determinant;

    // First row of inverse matrix
    inverse[0] = (double) (matrix[3] * matrix[5] - (matrix[4] * matrix[4]));
    inverse[1] = (double) -(matrix[1] * matrix[5] - matrix[4] * matrix[2]);
    inverse[2] = (double) (matrix[1] * matrix[4] - matrix[3] * matrix[2]);

    // Calculate determinant.
    determinant = (double) matrix[0] * inverse[0] + matrix[1] * inverse[1] + matrix[2] * inverse[2];

    // Apply determinant of the first row.
    inverse[0] /= determinant;
    inverse[1] /= determinant;
    inverse[2] /= determinant;

    // Second row of inverse matrix
    inverse[3] = ((double) (matrix[0] * matrix[5] - (matrix[2] * matrix[2]))) / determinant;
    inverse[4] = ((double) -(matrix[0] * matrix[4] - matrix[1] * matrix[2])) / determinant;

    // Third row of inverse matrix
    inverse[5] = ((double) (matrix[0] * matrix[3] - (matrix[1] * matrix[1]))) / determinant;
}

int64_t matrixInverseSymmetric3by3_int(int32_t *matrix, int64_t *inverse)
{
    // First row of adjoint matrix
    inverse[0] = (int64_t) (matrix[3] * matrix[5] - (matrix[4] * matrix[4]));
    inverse[1] = (int64_t) -(matrix[1] * matrix[5] - matrix[4] * matrix[2]);
    inverse[2] = (int64_t) (matrix[1] * matrix[4] - matrix[3] * matrix[2]);

    // Second row of adjoint matrix
    inverse[3] = (int64_t) (matrix[0] * matrix[5] - (matrix[2] * matrix[2]));
    inverse[4] = (int64_t) -(matrix[0] * matrix[4] - matrix[1] * matrix[2]);

    // Third row of adjoint matrix
    inverse[5] = (int64_t) (matrix[0] * matrix[3] - (matrix[1] * matrix[1]));

    // Return determinant.
    return (matrix[0] * inverse[0] + matrix[1] * inverse[1] + matrix[2] * inverse[2]);
}

void calculateTMatrix_double(int32_t *sm, double *s3i, double *T)
{
	// First row of T matrix
	T[0] = -(s3i[0] * sm[3] + s3i[1] * sm[4] + s3i[2] * sm[5]);
	T[1] = -(s3i[0] * sm[4] + s3i[1] * sm[9] + s3i[2] * sm[10]);
	T[2] = -(s3i[0] * sm[9] + s3i[1] * sm[13] + s3i[2] * sm[14]);

	// Second row of T matrix
	T[3] = -(s3i[1] * sm[3] + s3i[3] * sm[4] + s3i[4] * sm[5]);
	T[4] = -(s3i[1] * sm[4] + s3i[3] * sm[9] + s3i[4] * sm[10]);
	T[5] = -(s3i[1] * sm[9] + s3i[3] * sm[13] + s3i[4] * sm[14]);

	// Third row of T matrix
	T[6] = -(s3i[2] * sm[3] + s3i[4] * sm[4] + s3i[5] * sm[5]);
	T[7] = -(s3i[2] * sm[4] + s3i[4] * sm[9] + s3i[5] * sm[10]);
	T[8] = -(s3i[2] * sm[9] + s3i[4] * sm[13] + s3i[5] * sm[14]);
}

void calculateTMatrix_int(int32_t *sm, int64_t *s3i, int64_t *T)
{
	// First row of T matrix
	T[0] = -(s3i[0] * sm[3] + s3i[1] * sm[4] + s3i[2] * sm[5]);
	T[1] = -(s3i[0] * sm[4] + s3i[1] * sm[9] + s3i[2] * sm[10]);
	T[2] = -(s3i[0] * sm[9] + s3i[1] * sm[13] + s3i[2] * sm[14]);

	// Second row of T matrix
	T[3] = -(s3i[1] * sm[3] + s3i[3] * sm[4] + s3i[4] * sm[5]);
	T[4] = -(s3i[1] * sm[4] + s3i[3] * sm[9] + s3i[4] * sm[10]);
	T[5] = -(s3i[1] * sm[9] + s3i[3] * sm[13] + s3i[4] * sm[14]);

	// Third row of T matrix
	T[6] = -(s3i[2] * sm[3] + s3i[4] * sm[4] + s3i[5] * sm[5]);
	T[7] = -(s3i[2] * sm[4] + s3i[4] * sm[9] + s3i[5] * sm[10]);
	T[8] = -(s3i[2] * sm[9] + s3i[4] * sm[13] + s3i[5] * sm[14]);
}

void calculateSecondTermOfM_double(int32_t *sm, double *T, double *M)
{
	// First row of second term of M matrix
	M[0] = (T[0] * sm[3] + T[3] * sm[4] + T[6] * sm[5]);
	M[1] = (T[1] * sm[3] + T[4] * sm[4] + T[7] * sm[5]);
	M[2] = (T[2] * sm[3] + T[5] * sm[4] + T[8] * sm[5]);

	// Second row of second term of M matrix
	M[3] = (T[0] * sm[4] + T[3] * sm[9] + T[6] * sm[10]);
	M[4] = (T[1] * sm[4] + T[4] * sm[9] + T[7] * sm[10]);
	M[5] = (T[2] * sm[4] + T[5] * sm[9] + T[8] * sm[10]);

	// Third row of second term of M matrix
	M[6] = (T[0] * sm[9] + T[3] * sm[13] + T[6] * sm[14]);
	M[7] = (T[1] * sm[9] + T[4] * sm[13] + T[7] * sm[14]);
	M[8] = (T[2] * sm[9] + T[5] * sm[13] + T[8] * sm[14]);
}

void calculateSecondTermOfM_int(int32_t *sm, int64_t *T, int32_t *M, int64_t determinant)
{
	// First row of second term of M matrix
	M[0] = (int32_t) ((T[0] * sm[3] + T[3] * sm[4] + T[6] * sm[5]) / determinant);
	M[1] = (int32_t) ((T[1] * sm[3] + T[4] * sm[4] + T[7] * sm[5]) / determinant);
	M[2] = (int32_t) ((T[2] * sm[3] + T[5] * sm[4] + T[8] * sm[5]) / determinant);

	// Second row of second term of M matrix
	M[3] = (int32_t) ((T[0] * sm[4] + T[3] * sm[9] + T[6] * sm[10]) / determinant);
	M[4] = (int32_t) ((T[1] * sm[4] + T[4] * sm[9] + T[7] * sm[10]) / determinant);
	M[5] = (int32_t) ((T[2] * sm[4] + T[5] * sm[9] + T[8] * sm[10]) / determinant);

	// Third row of second term of M matrix
	M[6] = (int32_t) ((T[0] * sm[9] + T[3] * sm[13] + T[6] * sm[14]) / determinant);
	M[7] = (int32_t) ((T[1] * sm[9] + T[4] * sm[13] + T[7] * sm[14]) / determinant);
	M[8] = (int32_t) ((T[2] * sm[9] + T[5] * sm[13] + T[8] * sm[14]) / determinant);
}

void calculateCompleteM_double(int32_t *sm, double *M)
{
	// Add the two matrices S1 and -S2*S3^{-1}*S2^{T} and perform
	// left multiplication with C^{-1} simultaneously.
	M[0] = (M[0] + sm[0]) / 2;
	M[1] = (M[1] + sm[1]) / 2;
	M[3] -= sm[1];
	M[2] = (M[2] + sm[2]) / 2;
	M[6] = (M[6] + sm[2]) / 2;
	M[4] -= sm[6];
	M[5] -= sm[7];
	M[7] = (M[7] + sm[7]) / 2;
	M[8] = (M[8] + sm[11]) / 2;
}

void calculateCompleteM_int(int32_t *sm, int32_t *M)
{
	// Add the two matrices S1 and -S2*S3^{-1}*S2^{T} and perform
	// left multiplication with C^{-1} simultaneously.
	M[0] = (M[0] + sm[0]) / 2;
	M[1] = (M[1] + sm[1]) / 2;
	M[3] -= sm[1];
	M[2] = (M[2] + sm[2]) / 2;
	M[6] = (M[6] + sm[2]) / 2;
	M[4] -= sm[6];
	M[5] -= sm[7];
	M[7] = (M[7] + sm[7]) / 2;
	M[8] = (M[8] + sm[11]) / 2;
}

void conicToGeneralForm_double(double *cc, EllipseDouble *results)
{
	static double centerPointDenominator = 0;
	static double mu = 0;
	// Storing to avoid recalculation.
	static double muA = 0;
	static double muB = 0;
	static double muC = 0;
	static double first_term = 0;
	static double sqrt_term = 0;

	centerPointDenominator = 2 * ((cc[1] * cc[1]) - (cc[0] * cc[2]));
	results->centerX = ((cc[2] * cc[3]) - (cc[1] * cc[4])) / centerPointDenominator;
	results->centerY = ((cc[0] * cc[3]) - (cc[1] * cc[3])) / centerPointDenominator;

	mu = 1 / ((cc[0] * results->centerX * results->centerX) +
		      (2 * cc[1] * results->centerX * results->centerY) +
			  (cc[2] * results->centerY * results->centerY) -
			  cc[5]);

	muA = mu * cc[0];
	muB = mu * cc[1];
	muC = mu * cc[2];

	sqrt_term = sqrt(((muA - muC) * (muA - muC)) + (4 * muB * muB)) / 2;
	first_term = (muA + muC) / 2;

	results->minorAxis = 1 / sqrt(first_term + sqrt_term);
	results->majorAxis = 1 / sqrt(first_term - sqrt_term);

	results->rotationAngle = atan2(-2 * cc[1], cc[2] - cc[0]);
	results->isValid = true;
}
