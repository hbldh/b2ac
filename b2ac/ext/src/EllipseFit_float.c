#include "EllipseFit_float.h"

EllipseFloat FitEllipse_float(int32_t *x, int32_t *y, int16_t numberOfPoints)
{
    int16_t _i;
    // Storage for the scatter matrix. It is a 6x6 symmetrical matrix.
    int64_t scatterMatrix[21] = {0};
    float T[9] = {0.0};
	// Eigensystem matrix storage.
	float M[9] = {0.0};
	// Eigenvectors and eigenvalues storage.
	float eigenvalues[3] =  {0.0};
	float eigenvector[3] =  {0.0};

	float conicCoefficients[6] = {0};

    int32_t _xMean = 0;
    int32_t _yMean = 0;
	int16_t ellipticalEigenvectorIndex = -1;
	int16_t eigenvalueOrderToTry[3] = {1, 0, 2};

	EllipseFloat results;

    _xMean = findMeanOfCoordinates(x, numberOfPoints);
    _yMean = findMeanOfCoordinates(y, numberOfPoints);

    // Calculate the T and M matrices needed.
    calculateScatterMatrix(x, y, _xMean, _yMean,
    		numberOfPoints, scatterMatrix);
    calculateMatricesMAndT_float(scatterMatrix, T, M);

	// Now we have M, apply eigensolver to get a1 vector.
	QRAlgorithm_float(M, eigenvalues);

	// Eigenvalues calculated; evaluate which one is elliptical.
	// It can be proven that only one eigenvalue will be positive in case
	// of elliptical solution; if none fulfills this, then fitting
	// was unsuccessful.

	for(_i = 0; _i < 3; _i++)
	{
		// Find the eigenvector to this eigenvalue.
		inverseIteration_float(M, eigenvector,
				eigenvalues[eigenvalueOrderToTry[_i]], 1);

		if (((4 * eigenvector[0] * eigenvector[2]) -
		    (eigenvector[1] * eigenvector[1])) > 0.0)
		{
			ellipticalEigenvectorIndex = eigenvalueOrderToTry[_i];
			break;
		}
		else
		{
			ellipticalEigenvectorIndex = -1;
		}
	}

	if (ellipticalEigenvectorIndex == -1)
	{
		// Non-elliptical solution found. Returning empty results.
		results.isValid = false;
	}
	else
	{
		// Extract the coefficients calculated and find the three last ones.
		// Flip the eigenvector to a positive direction in case a negative
		// one has been found. This makes the rotation of the ellipse conform.
		if (eigenvector[0] < 0)
		{
			conicCoefficients[0] = -eigenvector[0];
			conicCoefficients[1] = -eigenvector[1];
			conicCoefficients[2] = -eigenvector[2];
		}
		else
		{
			conicCoefficients[0] = eigenvector[0];
			conicCoefficients[1] = eigenvector[1];
			conicCoefficients[2] = eigenvector[2];
		}
		conicCoefficients[3] = conicCoefficients[0] * T[0] +
							   conicCoefficients[1] * T[1] +
							   conicCoefficients[2] * T[2];
		conicCoefficients[4] = conicCoefficients[0] * T[3] +
							   conicCoefficients[1] * T[4] +
							   conicCoefficients[2] * T[5];
		conicCoefficients[5] = conicCoefficients[0] * T[6] +
							   conicCoefficients[1] * T[7] +
							   conicCoefficients[2] * T[8];

		// Convert to general form and store this information
		// in the output struct.
		conicToGeneralForm_float(conicCoefficients, &results);

		// Add the previously subtracted mean point.
		results.centerX += _xMean;
		results.centerY += _yMean;
	}
	return (results);
}

void calculateMatricesMAndT_float(int64_t *scatterMatrix, float *T, float *M)
{
    // Storage for matrices to calculate on the way to eigensystem matrix.
    float S3Inverse[6] = {0.0};

    matrixInverseSymmetric3by3_float(&scatterMatrix[15], S3Inverse);
    calculateTMatrix_float(scatterMatrix, S3Inverse, T);
	calculateSecondTermOfM_float(scatterMatrix, T, M);
	calculateCompleteM_float(scatterMatrix, M);
}

void matrixInverseSymmetric3by3_float(int64_t *matrix, float *inverse)
{
    float determinant;

    // First row of inverse matrix
    inverse[0] = (float) (matrix[3] * matrix[5] - (matrix[4] * matrix[4]));
    inverse[1] = (float) -(matrix[1] * matrix[5] - matrix[4] * matrix[2]);
    inverse[2] = (float) (matrix[1] * matrix[4] - matrix[3] * matrix[2]);

    // Calculate determinant.
    determinant = (float) matrix[0] * inverse[0] + matrix[1] * inverse[1] + matrix[2] * inverse[2];

    // Apply determinant of the first row.
    inverse[0] /= determinant;
    inverse[1] /= determinant;
    inverse[2] /= determinant;

    // Second row of inverse matrix
    inverse[3] = ((float) (matrix[0] * matrix[5] - (matrix[2] * matrix[2]))) / determinant;
    inverse[4] = ((float) -(matrix[0] * matrix[4] - matrix[1] * matrix[2])) / determinant;

    // Third row of inverse matrix
    inverse[5] = ((float) (matrix[0] * matrix[3] - (matrix[1] * matrix[1]))) / determinant;
}

void calculateTMatrix_float(int64_t *sm, float *s3i, float *T)
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

void calculateSecondTermOfM_float(int64_t *sm, float *T, float *M)
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

void calculateCompleteM_float(int64_t *sm, float *M)
{
	// Add the two matrices S1 and -S2*S3^{-1}*S2^{T} and perform
	// left multiplication with C^{-1} simultaneously.
	M[0] = (M[0] + sm[0]) / 2;
	M[1] = (M[1] + sm[1]) / 2;
	M[3] = -(M[3] + sm[1]);
	M[2] = (M[2] + sm[2]) / 2;
	M[6] = (M[6] + sm[2]) / 2;
	M[4] = -(M[4] + sm[6]);
	M[5] = -(M[5] + sm[7]);
	M[7] = (M[7] + sm[7]) / 2;
	M[8] = (M[8] + sm[11]) / 2;
}

void conicToGeneralForm_float(float *cc,
		EllipseFloat *results)
{
    results->rotationAngle = atan2f(cc[1], cc[0] - cc[2]) / 2;

    float cost = cosf(results->rotationAngle);
    float sint = sinf(results->rotationAngle);

	// Using indirect trigonometric solution.
    /*
    float givensValues[2] = {0.0};
    givensRotation_float(givensValues, cc[1], cc[0] - cc[2]);

    if (givensValues[1] > 0)
    {
        sint = sqrtf((1 - (-givensValues[1]))/2);
        cost = givensValues[0] / (2 * sint);
    }
    else
    {
        cost = sqrtf(((-givensValues[1]) + 1)/2);
        sint = givensValues[0] / (2 * cost);
    }
    */

    float cost_sq = cost * cost;
    float sint_sq = sint * sint;

    float a_prime = (cc[0] * cost_sq) +
    				 (cc[1] * cost * sint) +
    				 (cc[2] * sint_sq);
    float c_prime = (cc[0] * sint_sq) -
    				 (cc[1] * cost * sint) +
					 (cc[2] * cost_sq);
    float d_prime = (cc[3] * cost) + (cc[4] * sint);
    float e_prime = (-(cc[3] * sint)) + (cc[4] * cost);

    float x_prime = (-d_prime) / (2 * a_prime);
    float y_prime = (-e_prime) / (2 * c_prime);

    results->centerX = (x_prime * cost) - (y_prime * sint);
    results->centerY = (x_prime * sint) + (y_prime * cost);

    float numerator = (-4 * cc[5] * a_prime * c_prime) +
                       (c_prime * (d_prime * d_prime)) +
                       (a_prime * (e_prime * e_prime));

    if (a_prime > c_prime)
    {
    	results->yAxis = sqrtf(numerator /
    		(4 * a_prime * (c_prime * c_prime)));
    	results->xAxis = sqrtf(numerator /
    		(4 * (a_prime * a_prime) * c_prime));
    }
    else
    {
    	results->xAxis = sqrtf(numerator /
    	    (4 * a_prime * (c_prime * c_prime)));
    	results->yAxis = sqrtf(numerator /
    	    (4 * (a_prime * a_prime) * c_prime));
    }

    results->isValid = true;
}
