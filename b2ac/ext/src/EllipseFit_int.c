#include "EllipseFit_int.h"

EllipseFloat FitEllipse_int(int32_t *x, int32_t *y, int16_t numberOfPoints)
{
    int16_t _i;
    // Storage for the scatter matrix. It is a 6x6 symmetrical matrix.
    int64_t scatterMatrix[21] = {0};
    int64_t S3Determinant = 0;
	int64_t T[9] = {0};
	// Eigensystem matrix storage.
	int64_t M[9] = {0};
	// Eigenvectors and eigenvalues storage.
	int64_t eigenvalues[3] =  {0};
	int64_t eigenvector[3] =  {0};
	static int64_t conicCoefficients[6] = {0};

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
    S3Determinant = calculateMatricesMAndT_int(scatterMatrix, T, M);

    // Now we have M, apply eigensolver to get a1 vector.
    QRAlgorithm_int(M, eigenvalues);

	// Eigenvalues calculated; evaluate which one is elliptical.
	// It can be proven that only one eigenvalue will be positive in case
	// of elliptical solution; if none fulfills this, then fitting was unsucessful.

	for(_i = 0; _i < 3; _i++)
	{
		// Find the eigenvector to this eigenvalue.
		inverseIteration_int(M, eigenvector,
				eigenvalues[eigenvalueOrderToTry[_i]]);

		if ((4 * eigenvector[0] * eigenvector[2]) -
				(eigenvector[1] * eigenvector[1]) > 0.0)
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
		conicCoefficients[3] = (conicCoefficients[0] * T[0] +
								conicCoefficients[1] * T[1] +
								conicCoefficients[2] * T[2]) / S3Determinant;
		conicCoefficients[4] = (conicCoefficients[0] * T[3] +
								conicCoefficients[1] * T[4] +
								conicCoefficients[2] * T[5]) / S3Determinant;
		conicCoefficients[5] = (conicCoefficients[0] * T[6] +
								conicCoefficients[1] * T[7] +
								conicCoefficients[2] * T[8]) / S3Determinant;

		// Convert to general form and store this information
		// in the output struct.
		conicToGeneralForm_int(conicCoefficients, &results);

		// Add the previously subtracted mean point.
		results.centerX += _xMean;
		results.centerY += _yMean;
	}

	return (results);
}

int64_t calculateMatricesMAndT_int(int64_t *scatterMatrix, int64_t *T, int64_t *M)
{
    int64_t _S3Determinant = 0;
    // Storage for variables and matrices to calculate on the way to eigensystem matrix.
    int64_t _S3Inverse[6] = {0};

    _S3Determinant = matrixInverseSymmetric3by3_int(&scatterMatrix[15], _S3Inverse);
    calculateTMatrix_int(scatterMatrix, _S3Inverse, T);
    _S3Determinant = scaleTMatrix(T, _S3Determinant);
	calculateSecondTermOfM_int(scatterMatrix, T, M, _S3Determinant);
	calculateCompleteM_int(scatterMatrix, M);

    return (_S3Determinant);
}

int64_t matrixInverseSymmetric3by3_int(int64_t *A, int64_t *inverse)
{
    // First row of adjoint matrix
    inverse[0] = (A[3] * A[5] - (A[4] * A[4]));
    inverse[1] = -(A[1] * A[5] - A[4] * A[2]);
    inverse[2] = (A[1] * A[4] - A[3] * A[2]);

    // Second row of adjoint matrix
    inverse[3] = (A[0] * A[5] - (A[2] * A[2]));
    inverse[4] = -(A[0] * A[4] - A[1] * A[2]);

    // Third row of adjoint matrix
    inverse[5] = (A[0] * A[3] - (A[1] * A[1]));

    // Return determinant.
    return (A[0] * inverse[0] + A[1] * inverse[1] + A[2] * inverse[2]);
}

void calculateTMatrix_int(int64_t *scatterMatrix, int64_t *s3Inverse,
		int64_t *T)
{
	// First row of T matrix
	T[0] = -(s3Inverse[0] * scatterMatrix[3] +
		 	 s3Inverse[1] * scatterMatrix[4] +
			 s3Inverse[2] * scatterMatrix[5]);
	T[1] = -(s3Inverse[0] * scatterMatrix[4] +
			 s3Inverse[1] * scatterMatrix[9] +
			 s3Inverse[2] * scatterMatrix[10]);
	T[2] = -(s3Inverse[0] * scatterMatrix[9] +
			 s3Inverse[1] * scatterMatrix[13] +
			 s3Inverse[2] * scatterMatrix[14]);

	// Second row of T matrix
	T[3] = -(s3Inverse[1] * scatterMatrix[3] +
			 s3Inverse[3] * scatterMatrix[4] +
			 s3Inverse[4] * scatterMatrix[5]);
	T[4] = -(s3Inverse[1] * scatterMatrix[4] +
			 s3Inverse[3] * scatterMatrix[9] +
			 s3Inverse[4] * scatterMatrix[10]);
	T[5] = -(s3Inverse[1] * scatterMatrix[9] +
			 s3Inverse[3] * scatterMatrix[13] +
			 s3Inverse[4] * scatterMatrix[14]);

	// Third row of T matrix
	T[6] = -(s3Inverse[2] * scatterMatrix[3] +
			 s3Inverse[4] * scatterMatrix[4] +
			 s3Inverse[5] * scatterMatrix[5]);
	T[7] = -(s3Inverse[2] * scatterMatrix[4] +
			 s3Inverse[4] * scatterMatrix[9] +
			 s3Inverse[5] * scatterMatrix[10]);
	T[8] = -(s3Inverse[2] * scatterMatrix[9] +
			 s3Inverse[4] * scatterMatrix[13] +
			 s3Inverse[5] * scatterMatrix[14]);
}

int64_t scaleTMatrix(int64_t *T, int64_t S3Determinant)
{
	int i;
	int64_t maxTValue = -1;
	int32_t scale = 0;

	for (i = 0; i < 9; i++)
	{
		if (labs(T[i]) > maxTValue)
		{
			maxTValue = labs(T[i]);
		}
	}
	// If maximal value is smaller than 2^32, don't do any scaling.
	if (maxTValue <= 4294967296L)
	{
		return (S3Determinant);
	}
	else if (maxTValue <= 1099511627776L)
	{
		// If maximal value is (2^32, 2^40], scale by 8 steps.
		scale = 8;
	}
	else if (maxTValue < 281474976710656L)
	{
		// If maximal value is (2^40, 2^48], scale by 16 steps.
		scale = 16;
	}
	else if (maxTValue < 72057594037927936L)
	{
		// If maximal value is (2^48, 2^56], scale by 24 steps.
		scale = 24;
	}
	else
	{
		// If maximal value is (2^56, 2^63], scale by 32 steps.
		scale = 32;
	}

	for (i = 0; i < 9; i++)
	{
		T[i] = T[i] >> scale;
	}
	return (S3Determinant >> scale);
}

void calculateSecondTermOfM_int(int64_t *scatterMatrix, int64_t *T,
		int64_t *M, int64_t determinant)
{
	// First row of second term of M matrix
	M[0] = ((T[0] * scatterMatrix[3] +
			 T[3] * scatterMatrix[4] +
			 T[6] * scatterMatrix[5]) / determinant);
	M[1] = ((T[1] * scatterMatrix[3] +
			 T[4] * scatterMatrix[4] +
			 T[7] * scatterMatrix[5]) / determinant);
	M[2] = ((T[2] * scatterMatrix[3] +
			 T[5] * scatterMatrix[4] +
			 T[8] * scatterMatrix[5]) / determinant);

	// Second row of second term of M matrix
	M[3] = ((T[0] * scatterMatrix[4] +
			 T[3] * scatterMatrix[9] +
			 T[6] * scatterMatrix[10]) / determinant);
	M[4] = ((T[1] * scatterMatrix[4] +
			 T[4] * scatterMatrix[9] +
			 T[7] * scatterMatrix[10]) / determinant);
	M[5] = ((T[2] * scatterMatrix[4] +
			 T[5] * scatterMatrix[9] +
			 T[8] * scatterMatrix[10]) / determinant);

	// Third row of second term of M matrix
	M[6] = ((T[0] * scatterMatrix[9] +
			 T[3] * scatterMatrix[13] +
			 T[6] * scatterMatrix[14]) / determinant);
	M[7] = ((T[1] * scatterMatrix[9] +
			 T[4] * scatterMatrix[13] +
			 T[7] * scatterMatrix[14]) / determinant);
	M[8] = ((T[2] * scatterMatrix[9] +
			 T[5] * scatterMatrix[13] +
			 T[8] * scatterMatrix[14]) / determinant);
}

void calculateCompleteM_int(int64_t *scatterMatrix, int64_t *M)
{
	// Add the two matrices S1 and -S2*S3^{-1}*S2^{T} and perform
	// left multiplication with C^{-1} simultaneously.
	M[0] = (M[0] + scatterMatrix[0]) / 2;
	M[1] = (M[1] + scatterMatrix[1]) / 2;
	M[3] = -(M[3] +scatterMatrix[1]);
	M[2] = (M[2] + scatterMatrix[2]) / 2;
	M[6] = (M[6] + scatterMatrix[2]) / 2;
	M[4] = -(M[4] + scatterMatrix[6]);
	M[5] = -(M[5] + scatterMatrix[7]);
	M[7] = (M[7] + scatterMatrix[7]) / 2;
	M[8] = (M[8] + scatterMatrix[11]) / 2;
}

void conicToGeneralForm_int(int64_t *cc, EllipseFloat *results)
{
	int64_t givensValues[3];
	int64_t sin_t, cos_t;
    // Has scaling unity.
    int64_t cos_t_squared, sin_t_squared, unity;
	int64_t squareComponentsDiff = cc[0] - cc[2];

	results->rotationAngle = (float) atan2f((float) cc[1], (float) (cc[0] - cc[2])) / 2;

	// First, check if any important vales are equal to zero; if so,
	// handle specially to avoid divide by zero errors and unneccessary
	// computations and scaling.
	if (cc[1] == 0)
	{
		if (squareComponentsDiff < 0)
		{
			cos_t = 0;
			sin_t = 1;
		}
		else
		{
			cos_t = 1;
			sin_t = 0;
		}
		cos_t_squared = 1;
		sin_t_squared = 0;
		unity = 1;
	}
	else if (squareComponentsDiff == 0)
	{
		if (cc[1] < 0)
		{
			cos_t = 5;
			sin_t = -5;
		}
		else
		{
			cos_t = 5;
			sin_t = 5;
		}
		cos_t_squared = 25;
		sin_t_squared = 25;
		unity = 50;
	}
	else
	{
		// Use the Givens rotation solution for obtaining the
		// double angle cos and sin values given the conic
		// coefficients.
		givensRotation_int(givensValues, cc[1], cc[0] - cc[2]);

		// Use the double angle formula that creates an actual addition inside
		// the square root call.
		if (givensValues[1] > 0)
		{
			sin_t = sqrt_int64((uint64_t)
				((givensValues[2]  + givensValues[1]) >> 1), false);
			cos_t = givensValues[0] / (2 * sin_t);
		}
		else
		{
			cos_t = sqrt_int64((uint64_t)
				(((-givensValues[1]) + givensValues[2]) >> 1), false);
			sin_t = givensValues[0] / (2 * cos_t);
		}
		// Now that we have cos and sin values we establish a new scaling with
        // unity value obtained through the sin^2(x) + cos^2(x) = 1
		// trigonometric identity.
		cos_t_squared = cos_t * cos_t;
		sin_t_squared = sin_t * sin_t;
		unity = cos_t_squared + sin_t_squared;
	}

	// Now rotate the conic coefficient representation of the ellipse to
	// the canonical representation, i.e. without cross-term B.
    // Ends up with scaling unity.
    int64_t a_prime = (cc[0] * cos_t_squared) +
					  (cc[1] * cos_t * sin_t) +
					  (cc[2] * sin_t_squared);
    int64_t c_prime = (cc[0] * sin_t_squared) -
					  (cc[1] * cos_t * sin_t) +
					  (cc[2] * cos_t_squared);
    // Ends up with scaling sqrt(unity)
    int64_t d_prime = (cc[3] * cos_t) + (cc[4] * sin_t);
    int64_t e_prime = (-(cc[3] * sin_t)) + (cc[4] * cos_t);

    // Now calculate center points in the rotated ellipse coordinate system.
    int64_t x_prime_num = (-d_prime);
    int64_t x_prime_denom = (2 * a_prime);

    int64_t y_prime_num = (-e_prime);
    int64_t y_prime_denom = (2 * c_prime);

    // Rotate the center points to original coordinate system.
	results->centerX = (((float)(x_prime_num * cos_t) / x_prime_denom) -
			((float)(y_prime_num * sin_t) / y_prime_denom));
	results->centerY = (((float)(x_prime_num * sin_t) / x_prime_denom) +
			((float)(y_prime_num * cos_t) / y_prime_denom));

	// Now we have to rescale the rotated conic coefficients down a bit,
	// since there will be a cubic factors in the axis calculations.
	int64_t sqrt_unity = sqrt_int64(unity, false);
	a_prime /= unity;
	c_prime /= unity;
	d_prime /= sqrt_unity;
	e_prime /= sqrt_unity;

	int64_t numerator = ((-4 * cc[5] * a_prime * c_prime) +
                        (c_prime * (d_prime * d_prime)) +
                        (a_prime * (e_prime * e_prime)));
	uint64_t tmp_val;

	if (a_prime > c_prime)
	{
		tmp_val = sqrt_int64((numerator /
	    		(4 * a_prime * (c_prime * c_prime))), true);
		results->yAxis = convertSqrtToFloat(tmp_val);

    	tmp_val = sqrt_int64((numerator /
    		(4 * (a_prime * a_prime) * c_prime)), true);
    	results->xAxis = convertSqrtToFloat(tmp_val);
	}
	else
	{
		tmp_val = sqrt_int64((numerator /
				(4 * a_prime * (c_prime * c_prime))), true);
		results->xAxis = convertSqrtToFloat(tmp_val);

		tmp_val = sqrt_int64((numerator /
			(4 * (a_prime * a_prime) * c_prime)), true);
		results->yAxis = convertSqrtToFloat(tmp_val);
	}
    results->isValid = true;
}

/**
 * Converts integer square root function output to a floating
 * point value.
 */
float convertSqrtToFloat(uint64_t val)
{
	return (((float)(val >> 32)) +
			((float) ((val & 0xFFFFFFFF))) / 0xFFFFFFFF);
}

