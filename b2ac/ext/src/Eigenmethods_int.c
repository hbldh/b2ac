/**
 * @file	Eigenmethods_int.cpp
 * @author	hbh
 * @date	Jul 25, 2013
 * @version	0.1
 *
 * \brief TBD.
 */

#include "Eigenmethods_int.h"

void QRAlgorithm_int(int64_t *A, int64_t *eigenvalues)
{
	int64_t mu;
	int k;
	int m = 2;
	int mMapping[3] = {0, 4, 8};
	int32_t scale = 0;

	int64_t B[9];
	// Copy matrix, since we destroy it in the process of
	// finding eigenvalues.
	for (k = 0; k < 9; k++)
	{
		B[k] = A[k];
	}

	scale = scaleMatrix(B);

	matrixToHessenbergForm_Givens_int(B);

	while (m > 0)
	{
		// Apply the shift to improve algorithm convergence.
		// The index (m * 2) refer to the diagonal entries of the
		// A. Set the shift value to the "lower right corner" value.
		mu = B[mMapping[m]];
		for (k = 0; k <= m; k++)
		{
			B[mMapping[k]] -= mu;
		}

		// Run the Givens rotation producing RQ from A.
		performQRStep_Givens_int(B, m);

		// Add the shift to the diagonal elements again.
		for (k = 0; k <= m; k++)
		{
			B[mMapping[k]] += mu;
		}

		// If A[m, m-1] is equal to zero, an eigenvalue can be
		// deemed to have been found. Scale it up according to the
		// scale we reduced the matrix by previously.
		if (labs(B[mMapping[m] - 1]) <= 1)
		{
			eigenvalues[m] = B[mMapping[m]] << scale;
			m--;
		}
	}
	// Two eigenvalues have been found in the loop above. The final one is
	// left in the A A top left corner. Scale it up according to the
	// scale we reduced the matrix by previously.
	eigenvalues[0] = B[0] << scale;

	// Also sort the eigenvalues with simple Bubblesort.
	// Use mu as temporary variable for sorting.
	if (eigenvalues[0] > eigenvalues[1])
	{
		mu = eigenvalues[1];
		eigenvalues[1] = eigenvalues[0];
		eigenvalues[0] = mu;
	}
	if (eigenvalues[1] > eigenvalues[2])
	{
		mu = eigenvalues[2];
		eigenvalues[2] = eigenvalues[1];
		eigenvalues[1] = mu;
	}
	if (eigenvalues[0] > eigenvalues[1])
	{
		mu = eigenvalues[1];
		eigenvalues[1] = eigenvalues[0];
		eigenvalues[0] = mu;
	}
}

/**
 * Converts a 3x3 general matrix to a Hessenberg form, using
 * Givens rotations.
 *
 * @brief In Hessenberg form, a matrix has the same eigenvalues
 * as in its original form, but it is upper triangular, with
 * non-zero values from the first subdiagonal and up.
 *
 */
void matrixToHessenbergForm_Givens_int(int64_t *A)
{
	int64_t givensValues[3] = {0, 0, 0};
	int i;
	int64_t tau_one = 0;
	int64_t tau_two = 0;

	// For a 3x3 A, we need only rotate the A[2, 2] = A[6] value
	// into the A. Thus we first do a rotation with rows [1,2] and the
	// with columns [1,2].
	givensRotation_int(givensValues, A[3], A[6]);

	// Perform the row rotations.
	for(i = 0; i < 3; i++)
	{
		tau_one = A[3 + i];
		tau_two = A[6 + i];
		A[3 + i] = (tau_one * givensValues[0] -
					tau_two * givensValues[1]) / givensValues[2];
		A[6 + i] = (tau_one * givensValues[1] +
					tau_two * givensValues[0]) / givensValues[2];
	}

	// Perform the column rotations.
	for(i = 0; i < 9; i += 3)
	{
		tau_one = A[i + 1];
		tau_two = A[i + 2];
		A[i + 1] = (tau_one * givensValues[0] -
					tau_two * givensValues[1]) / givensValues[2];
		A[i + 2] = (tau_one * givensValues[1] +
					tau_two * givensValues[0]) / givensValues[2];
	}
}

/**
 * Performs QR algorithm iteration using Givens rotation.
 *
 * @brief The QR algorithm's key step is to perform a QR decomposition
 * of a Hessenberg matrix and multiply the matrices in reverse order,
 * i.e. RQ. The Givens solution here allows for implicit decomposition
 * and multiplication, yielding a much faster solution.
 *
 */
void performQRStep_Givens_int(int64_t *A, int m)
{
	int k, i;
	int row_ind, giv_ind;
	int64_t tau_one = 0;
	int64_t tau_two = 0;
	int64_t givensValues[6] = {0, 0, 0, 0, 0, 0};

	// If m == 2 then we are still working with the 3x3 matrix.
	if (m == 2)
	{
		// Perform row rotation operations, first with the [0,1] rows,
		// then with [1,2].
		for (k = 0; k < 2; k++)
		{
			// Set indices for obtaining correct data for this iteration.
			giv_ind = 3 * k;
			row_ind = 3 * k;
			// Find the Givens rotation to apply here.
			givensRotation_int(&givensValues[giv_ind],
									 A[row_ind + k],
									 A[row_ind + 3 + k]);
			for (i = 0; i < 3; i++)
			{
				tau_one = A[row_ind + i];
				tau_two = A[row_ind + 3 + i];
				A[row_ind + i] = (tau_one * givensValues[giv_ind] -
								  tau_two * givensValues[giv_ind + 1]) /
								  givensValues[giv_ind + 2];
				A[row_ind + 3 + i] = (tau_one * givensValues[giv_ind + 1] +
									  tau_two * givensValues[giv_ind]) /
									  givensValues[giv_ind + 2];
			}
		}
		// Perform column rotation operations.
		for (k = 0; k < 2; k++)
		{
			// Set indices for obtaining correct data for this iteration.
			giv_ind = 3 * k;
			for (i = 0; i < 3; i++)
			{
				row_ind = 3 * i;
				tau_one = A[row_ind + k];
				tau_two = A[row_ind + k + 1];
				A[row_ind + k] = (tau_one * givensValues[giv_ind] -
								  tau_two * givensValues[giv_ind + 1]) /
								  givensValues[giv_ind + 2];
				A[row_ind + k + 1] = (tau_one * givensValues[giv_ind + 1] +
									  tau_two * givensValues[giv_ind]) /
		  	  	  	  	  	  	  	  givensValues[giv_ind + 2];
			}
		}
	}
	else // m == 1, we will be using the upper left 2x2 A.
	{
		// Find the Givens rotation to apply here.
		givensRotation_int(givensValues, A[0], A[3]);

		// Perform row rotation operations.
		tau_one = A[0];
		tau_two = A[3];
		A[0] = (tau_one * givensValues[0] -
				tau_two * givensValues[1]) / givensValues[2];
		A[3] = (tau_one * givensValues[1] +
				tau_two * givensValues[0]) / givensValues[2];
		tau_one = A[1];
		tau_two = A[4];
		A[1] = (tau_one * givensValues[0] -
				tau_two * givensValues[1]) / givensValues[2];
		A[4] = (tau_one * givensValues[1] +
			    tau_two * givensValues[0]) / givensValues[2];

		// Perform column rotation operations.
		tau_one = A[0];
		tau_two = A[1];
		A[0] = (tau_one * givensValues[0] -
				tau_two * givensValues[1]) / givensValues[2];
		A[1] = (tau_one * givensValues[1] +
				tau_two * givensValues[0]) / givensValues[2];
		tau_one = A[3];
		tau_two = A[4];
		A[3] = (tau_one * givensValues[0] -
				tau_two * givensValues[1]) / givensValues[2];
		A[4] = (tau_one * givensValues[1] +
				tau_two * givensValues[0]) / givensValues[2];
	}
}

/**
 * Gets the cosine and sine values needed for performing Givens rotation.
 *
 * @brief Givens rotations are used for creating zeros in matrices
 * by rotating columns and/or rows of the matrix to accommodate for the
 * value in the cell that is zeroed by orthogonal matrix products.
 *
 */
void givensRotation_int(int64_t *csVect, int64_t a, int64_t b)
{
	csVect[0] = a;
	csVect[1] = -b;
	csVect[2] = sqrt_int64((a * a) + (b * b), false);
}

/**
 * 64 bit integer square root algorithm.
 *
 * @brief The calculations are the base-two analogue of the square
 * root algorithm we all learned in grammar school.  Since we're
 * in base 2, there is only one nontrivial trial multiplier.
 *
 * Notice that absolutely no multiplications or divisions are performed.
 * This means it'll be fast on a wide range of processors.
 *
 * Returns the value y s.t. max_{y} y * y <= x or, more programmatically,
 * floor(sqrt(x) * pow(2, 32))
 *
 */
int64_t sqrt_int64(uint64_t x, bool returnFractions) {
	int64_t a = 0L;                   // Accumulator
	int64_t r = 0L;                   // Remainder
	int64_t e = 0L;                   // Trial product

	int i;

	for (i = 0; i < 64; i++)
	{
		r = (r << 2) + ISQRT_TOP2BITS(x);
		x <<= 2;
		a <<= 1;
		e = (a << 1) + 1;
		if (r >= e)
		{
			  r -= e;
			  a++;
		}
	}

	if (returnFractions)
	{
		return (a);
	}
	else
	{
		return (a >> 32);
	}
}

/**
 * Inverse iteration for finding eigenvector.
 *
 * @brief Given a matrix and a known eigenvalue, this function
 * multiplies the matrix (A - mu*I)^{-1} with a iteratively improved
 * eigenvector approximation.
 *
 * Note that for a accurate eigenvalue, only one iteration is required for
 * a reasonable approximation.
 *
 * References:
 * [1] http://en.wikipedia.org/wiki/Inverse_iteration
 * [2] Gene H. Golub and Charles F. Van Loan. Matrix computations (3rd
 *     ed.). Johns Hopkins University Press, Baltimore, MD, USA, 1996.
 *
 */
void inverseIteration_int(int64_t *A, int64_t* eigenvector,
		int64_t eigenvalue)
{
	// Storage room for the inverse to form.
	int i;
	int64_t inverse[9];
	int64_t A_copy[9];

	// Make copy of matrix sent in, since we will subtract
	// values and then scale down, which is an irreversible operation.
	for (i = 0; i < 9; i++)
	{
		A_copy[i] = A[i];
	}

	// Subtract the eigenvalue from the matrix.
	A_copy[0] -= eigenvalue;
	A_copy[4] -= eigenvalue;
	A_copy[8] -= eigenvalue;

	// Downscale matrix and then invert it.
	scaleMatrix(A_copy);
	matrixInverse3by3_int(A_copy, inverse);

	// Perform matrix-vector multiplication.
	// In integer precision, use a [1, 1, 1]^T vector
	eigenvector[0] = inverse[0] + inverse[1] + inverse[2];
	eigenvector[1] = inverse[3] + inverse[4] + inverse[5];
	eigenvector[2] = inverse[6] + inverse[7] + inverse[8];

	// Scale eigenvector to reasonable size. We are only interested in the
	// relationship between the values not their actual values.
	// We rather use this special scaling than dividing by the actual
	// adjugate determinant, since we lose too much precision otherwise.
	// The returned scale is unimportant as well, since we only care about
	// the norm of the actual vector.
	scaleVector(eigenvector);
}

/**
 * Inverse of a general 3 by 3 matrix.
 */
int64_t matrixInverse3by3_int(int64_t *A, int64_t *inverse)
{
    // First row of adjugate A.
    inverse[0] = (A[4] * A[8] - A[7] * A[5]);  // Det #0
    inverse[1] = -(A[1] * A[8] - A[7] * A[2]);
    inverse[2] = (A[1] * A[5] - A[4] * A[2]);

    // Second row of adjugate A.
    inverse[3] = -(A[3] * A[8] - A[6] * A[5]);  // Det #1
    inverse[4] = (A[0] * A[8] - A[6] * A[2]);
    inverse[5] = -(A[0] * A[5] - A[3] * A[2]);

    // Third row of adjugate A.
    inverse[6] = (A[3] * A[7] - A[6] * A[4]);  // Det #2
    inverse[7] = -(A[0] * A[7] - A[6] * A[1]);
    inverse[8] = (A[0] * A[4] - A[3] * A[1]);

    // Return determinant.
    return (A[0] * inverse[0] + A[1] * inverse[3] + A[2] * inverse[6]);
}

/**
 * Scaling for QR algorithm matrices.
 *
 * @brief Simple, adaptive scaling for the matrix to find eigenvalues in.
 * It does yield more approximate eigenvalues, but for the ellipse fitting
 * algorithm, this does not matter very much.
 *
 */
int32_t scaleMatrix(int64_t *A)
{
	int i;
	int64_t maxTValue = -1;
	int32_t scale = 0;

	for (i = 0; i < 9; i++)
	{
		if (labs(A[i]) > maxTValue)
		{
			maxTValue = abs(A[i]);
		}
	}
	// If maximal value is smaller than 2^24, don't do any scaling.
	if (maxTValue <= 16777216L)
	{
		return (scale);
	}
	else if (maxTValue <= 4294967296L)
	{
		// If maximal value is (2^24, 2^32], scale by 8 steps.
		scale = 8;
	}
	else if (maxTValue <= 1099511627776L)
	{
		// If maximal value is (2^32, 2^40], scale by 16 steps.
		scale = 16;
	}
	else if (maxTValue < 281474976710656L)
	{
		// If maximal value is (2^40, 2^48], scale by 24 steps.
		scale = 24;
	}
	else if (maxTValue < 72057594037927936L)
	{
		// If maximal value is (2^48, 2^56], scale by 32 steps.
		scale = 32;
	}
	else
	{
		// If maximal value is (2^56, 2^63], scale by 40 steps.
		scale = 40;
	}

	for (i = 0; i < 9; i++)
	{
		A[i] = A[i] >> scale;
	}
	return (scale);
}

/**
 * Scaling for the eigenvectors generated by inverseIteration_int.
 *
 * @brief Simple, adaptive scaling for the eigenvector. This ensures
 * that all values inside the array is in the range [-2^16, 2^16].
 *
 *
 */
int32_t scaleVector(int64_t *v)
{
	int i;
	int64_t maxTValue = -1;
	int32_t scale = 0;

	for (i = 0; i < 3; i++)
	{
		if (labs(v[i]) > maxTValue)
		{
			maxTValue = labs(v[i]);
		}
	}

	// For every bit above 12, increase the scale by one to reduce the largest
	// value of the sent in vector to maximally +/- 4095.
	maxTValue = maxTValue >> 12;
	while(maxTValue > 0)
	{
	    maxTValue = maxTValue >> 1;
	    scale++;
	}

	for (i = 0; i < 3; i++)
	{
		v[i] = v[i] >> scale;
	}
	return (scale);
}
