#include "Eigenmethods_double.h"

/**
 * The QR Algorithm, for finding eigenvalues.
 *
 * @brief This version of the QR algorithms is implemented
 * purely using Givens rotations. It fills the input
 * vector eigenvalues with the eigenvalues it detects.
 */
void QRAlgorithm_double(double *A, double *eigenvalues)
{
	double mu;
	int k;
	int m = 2;
	int mMapping[3] = {0, 4, 8};

	double B[9];
	// Copy matrix, since we destroy it in the process of
	// finding eigenvalues.
	for (k = 0; k < 9; k++)
	{
		B[k] = A[k];
	}

	// First, convert the A to Hessenberg form.
	matrixToHessenbergForm_Givens_double(B);

	while (m > 0)
	{
		// Apply the shift to improve algorithm convergence.
		// The index (m * 2) refer to the diagonal entries of the
		// B. Set the shift value to the "lower right corner" value.
		mu = B[mMapping[m]];
		for (k = 0; k <= m; k++)
		{
			B[mMapping[k]] -= mu;
		}

		// Run the Givens rotation producing RQ from B.
		performQRStep_Givens_double(B, m);

		// Add the shift to the diagonal elements again.
		for (k = 0; k <= m; k++)
		{
			B[mMapping[k]] += mu;
		}

		// If B[m, m-1] is sufficiently small, an eigenvalue can be
		// deemed to have been found.

		if (fabs(B[mMapping[m] - 1]) < QR_ALGORITHM_TOLERANCE_DOUBLE)
		{
			eigenvalues[m] = B[mMapping[m]];
			m--;
		}
	}

	// Two eigenvalues have been found in the loop above. The final one is
	// left in the A A top left corner.
	eigenvalues[0] = B[0];

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
void matrixToHessenbergForm_Givens_double(double *A)
{
	double givensValues[2] = {0.0, 0.0};
	int i;
	double tau_one = 0.0;
	double tau_two = 0.0;

	// For a 3x3 A, we need only rotate the A[2, 2] = A[6] value
	// into the A. Thus we first do a rotation with rows [1,2] and the
	// with columns [1,2].
	givensRotation_double(givensValues, A[3], A[6]);

	// Perform the row rotations.
	for(i = 0; i < 3; i++)
	{
		tau_one = A[3 + i];
		tau_two = A[6 + i];
		A[3 + i] = tau_one * givensValues[0] - tau_two * givensValues[1];
		A[6 + i] = tau_one * givensValues[1] + tau_two * givensValues[0];
	}

	// Perform the column rotations.
	for(i = 0; i < 9; i += 3)
	{
		tau_one = A[i + 1];
		tau_two = A[i + 2];
		A[i + 1] = tau_one * givensValues[0] - tau_two * givensValues[1];
		A[i + 2] = tau_one * givensValues[1] + tau_two * givensValues[0];
	}
}

/**
 * Converts a 3x3 general matrix to a Hessenberg form, using
 * Householder reflections.
 *
 * @brief In Hessenberg form, a matrix has the same eigenvalues
 * as in its original form, but it is upper triangular, with
 * non-zero values from the first subdiagonal and up.
 *
 */
void matrixToHessenbergForm_Householder_double(double *A)
{
    double vk[2] = {0.0, 0.0,};
    double vk_norm = 0.0;
    double dot_array[3] = {0.0, 0.0, 0.0};

    // First iteration, fix zero at position (3,1).
    // Find Householder reflector.
    vk[0] = A[3];
    vk[1] = A[6];
    vk_norm = sqrt(vk[0] * vk[0] + vk[1] * vk[1]);
    vk[0] += signbit(vk[0]) * vk_norm;
    vk_norm = sqrt(vk[0] * vk[0] + vk[1] * vk[1]);
    vk[0] /= vk_norm;
    vk[1] /= vk_norm;

    // Apply Householder reflector.
    dot_array[0] = vk[0] * A[3] + vk[1] * A[6];
    dot_array[1] = vk[0] * A[4] + vk[1] * A[7];
    dot_array[2] = vk[0] * A[5] + vk[1] * A[8];
    A[3] -= 2 * vk[0] * dot_array[0];
    A[4] -= 2 * vk[0] * dot_array[1];
    A[5] -= 2 * vk[0] * dot_array[2];
    A[6] -= 2 * vk[1] * dot_array[0];
    A[7] -= 2 * vk[1] * dot_array[1];
    A[8] -= 2 * vk[1] * dot_array[2];

    dot_array[0] = vk[0] * A[1] + vk[1] * A[2];
    dot_array[1] = vk[0] * A[4] + vk[1] * A[5];
    dot_array[2] = vk[0] * A[7] + vk[1] * A[8];
    A[1] -= 2 * vk[0] * dot_array[0];
    A[2] -= 2 * vk[1] * dot_array[0];
    A[4] -= 2 * vk[0] * dot_array[1];
    A[5] -= 2 * vk[1] * dot_array[1];
    A[7] -= 2 * vk[0] * dot_array[2];
    A[8] -= 2 * vk[1] * dot_array[2];

    // Perform second and final iteration.
    // Find Householder reflector.
    A[2] = -A[2];
    A[5] = -A[5];
    A[7] = -A[7];

    // Done.
}

/**
 * Gets the cosine and sine values needed for performing Givens rotation.
 *
 * @brief Givens rotations are used for creating zeros in matrices
 * by rotating columns and/or rows of the matrix to accommodate for the
 * value in the cell that is zeroed by orthogonal matrix products.
 *
 */
void givensRotation_double(double *csVect, double a, double b)
{
	double tau;
	double tmpVal;
	if (b == 0.0)
	{
		csVect[0] = 1.0;
		csVect[1] = 0.0;
	}
	else
	{
		if (fabs(b) > fabs(a))
		{
			tau = (-a) / b;
			tmpVal = sqrt(1 + (tau * tau));
			csVect[1] = 1 / tmpVal;
			csVect[0] = tau / tmpVal;
		}
		else
		{
			tau = (-b) / a;
			tmpVal = sqrt(1 + (tau * tau));
			csVect[0] = 1 / tmpVal;
			csVect[1] = tau / tmpVal;
		}
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
void performQRStep_Givens_double(double *A, int m)
{
	int k, i;
	int row_ind, giv_ind;
	double tau_one = 0.0;
	double tau_two = 0.0;
	double givensValues[4] = {0.0, 0.0, 0.0, 0.0};

	// If m == 2 then we are still working with the 3x3 matrix.
	if (m == 2)
	{
		// Perform row rotation operations, first with the [0,1] rows,
		// then with [1,2].
		for (k = 0; k < 2; k++)
		{
			// Set indices for obtaining correct data for this iteration.
			giv_ind = 2 * k;
			row_ind = 3 * k;
			// Find the Givens rotation to apply here.
			givensRotation_double(&givensValues[giv_ind],
									 A[row_ind + k],
									 A[row_ind + 3 + k]);
			for (i = 0; i < 3; i++)
			{
				tau_one = A[row_ind + i];
				tau_two = A[row_ind + 3 + i];
				A[row_ind + i] = tau_one * givensValues[giv_ind] -
								 tau_two * givensValues[giv_ind + 1];
				A[row_ind + 3 + i] = tau_one * givensValues[giv_ind + 1] +
									 tau_two * givensValues[giv_ind];
			}
		}
		// Perform column rotation operations.
		for (k = 0; k < 2; k++)
		{
			// Set indices for obtaining correct data for this iteration.
			giv_ind = 2 * k;
			for (i = 0; i < 3; i++)
			{
				row_ind = 3 * i;
				tau_one = A[row_ind + k];
				tau_two = A[row_ind + k + 1];
				A[row_ind + k] = tau_one * givensValues[giv_ind] -
								 tau_two * givensValues[giv_ind + 1];
				A[row_ind + k + 1] = tau_one * givensValues[giv_ind + 1] +
									 tau_two * givensValues[giv_ind];
			}
		}
	}
	else // m == 1, we will be using the upper left 2x2 A.
	{
		// Find the Givens rotation to apply here.
		givensRotation_double(givensValues, A[0], A[3]);

		// Perform row rotation operations.
		tau_one = A[0];
		tau_two = A[3];
		A[0] = tau_one * givensValues[0] -
			   tau_two * givensValues[1];
		A[3] = tau_one * givensValues[1] +
			   tau_two * givensValues[0];
		tau_one = A[1];
		tau_two = A[4];
		A[1] = tau_one * givensValues[0] -
			   tau_two * givensValues[1];
		A[4] = tau_one * givensValues[1] +
			   tau_two * givensValues[0];

		// Perform column rotation operations.
		tau_one = A[0];
		tau_two = A[1];
		A[0] = tau_one * givensValues[0] -
			   tau_two * givensValues[1];
		A[1] = tau_one * givensValues[1] +
			   tau_two * givensValues[0];
		tau_one = A[3];
		tau_two = A[4];
		A[3] = tau_one * givensValues[0] -
			   tau_two * givensValues[1];
		A[4] = tau_one * givensValues[1] +
			   tau_two * givensValues[0];
	}
}

/**
 * Inverse iteration for finding eigenvector.
 *
 * @brief Given a matrix and a known eigenvalue, this function
 * multiplies the matrix (A - mu*I)^{-1} with a iteratively improved
 * eigenvector approximation.
 *
 * Note that for a accurate eigenvalue, one one iteration is required for
 * a reasonable approximation.
 *
 * References:
 * [1] http://en.wikipedia.org/wiki/Inverse_iteration
 * [2] Gene H. Golub and Charles F. Van Loan. Matrix computations (3rd
 *     ed.). Johns Hopkins University Press, Baltimore, MD, USA, 1996.
 *
 */
void inverseIteration_double(double *A, double* eigenvector,
		double eigenvalue, int32_t numberOfIterations)
{
	// Storage room for the inverse to form.
	int i;
	double determinant;
	double inverse[9];
	double tmp_storage[3];
	double norm_value = sqrt((double) 3);

	// Reset the eigenvector to a uniform unit vector.
	eigenvector[0] = 1 / norm_value;
	eigenvector[1] = eigenvector[0];
	eigenvector[2] = eigenvector[0];

	// Subtract the eigenvalue from the matrix.
	A[0] -= eigenvalue;
	A[4] -= eigenvalue;
	A[8] -= eigenvalue;

    // By not dividing by the determinant, but letting the normalization
    // in the end take care of that "scaling parameter", we are
    // more numerically stable against singular matrices.
	determinant = matrixInverse3by3_double(A, inverse);

	for (i = 0; i < numberOfIterations; i++)
	{
		// Perform matrix-vector multiplication.
		tmp_storage[0] = inverse[0] * eigenvector[0] +
						 inverse[1] * eigenvector[1] +
						 inverse[2] * eigenvector[2];
		tmp_storage[1] = inverse[3] * eigenvector[0] +
						 inverse[4] * eigenvector[1] +
						 inverse[5] * eigenvector[2];
		tmp_storage[2] = inverse[6] * eigenvector[0] +
						 inverse[7] * eigenvector[1] +
						 inverse[8] * eigenvector[2];

		// Normalise the eigenvector candidate.
		norm_value = sqrt(tmp_storage[0] * tmp_storage[0] +
						  tmp_storage[1] * tmp_storage[1] +
						  tmp_storage[2] * tmp_storage[2]);
		eigenvector[0] = tmp_storage[0] / norm_value;
		eigenvector[1] = tmp_storage[1] / norm_value;
		eigenvector[2] = tmp_storage[2] / norm_value;
	}

	// Add the eigenvalue back to the matrix.
	A[0] += eigenvalue;
	A[4] += eigenvalue;
	A[8] += eigenvalue;
}

/**
 * Inverse of a general 3 by 3 matrix.
 *
 * N.B. it returns the adjugate matrix and determinant separately
 * in order to be able to apply it only if it is needed.
 *
 */
double matrixInverse3by3_double(double *A, double *inverse)
{
    // First row of adjugate A
    inverse[0] = (A[4] * A[8] - A[7] * A[5]);  // Det #0
    inverse[1] = -(A[1] * A[8] - A[7] * A[2]);
    inverse[2] = (A[1] * A[5] - A[4] * A[2]);

    // Second row of adjugate A
    inverse[3] = -(A[3] * A[8] - A[6] * A[5]);  // Det #1
    inverse[4] = (A[0] * A[8] - A[6] * A[2]);
    inverse[5] = -(A[0] * A[5] - A[3] * A[2]);

    // Third row of adjugate A
    inverse[6] = (A[3] * A[7] - A[6] * A[4]);  // Det #2
    inverse[7] = -(A[0] * A[7] - A[6] * A[1]);
    inverse[8] = (A[0] * A[4] - A[3] * A[1]);

    return (A[0] * inverse[0] + A[1] * inverse[3] + A[2] * inverse[6]);
}
