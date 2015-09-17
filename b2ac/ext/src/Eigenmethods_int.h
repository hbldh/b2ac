/**
 * @file	Eigenmethods_int.h
 * @author	hbh
 * @date	Jul 25, 2013
 * @version	0.1
 *
 * \brief TBD.
 */

#ifndef EIGENMETHODS_INT_H_
#define EIGENMETHODS_INT_H_

#include <stdint.h>
#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
// Do not add anything in Windows.
#else
#include <stdbool.h>
#endif

// Help macro for the integer square root.
#define ISQRT_TOP2BITS(x) ((x & (3UL << (64 - 2))) >> (64 - 2))

void QRAlgorithm_int(int64_t *A, int64_t *eigenvalues);

void matrixToHessenbergForm_Givens_int(int64_t *A);

void performQRStep_Givens_int(int64_t *A, int m);

void givensRotation_int(int64_t *csVect, int64_t a, int64_t b);

int32_t scaleMatrix(int64_t *A);

int32_t scaleVector(int64_t *v);

int64_t sqrt_int64(uint64_t x, bool returnIntegerValue);

void inverseIteration_int(int64_t *A, int64_t* eigenvector,
		int64_t eigenvalue);

int64_t matrixInverse3by3_int(int64_t *A, int64_t *inverse);

#endif /* EIGENMETHODS_INT_H_ */
