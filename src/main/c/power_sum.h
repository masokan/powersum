#ifndef POWERSUM_H
#define POWERSUM_H

#include <gmp.h>
#include <stdio.h>
#include <time.h>

typedef struct PowerSum {
  /* To get the coefficients in the sum formula.
   * Parameters:
   *   power - desired power (IN)
   *   coeffList - an array coefficients (OUT)
   *   numCoeffs - number of elements in the array (OUT)
   * The returned array should be freed using freeCoefficients to avoid any
   * memory leak.
   */
  void (*getCoefficients)(long power, mpq_t ** coeffList, long * numCoeffs);

  /* To print the sum formula for a given power
   * Parameters:
   *   power - desired power (IN)
   *   out - file stream where the formula should be printed (IN)
   */
  void (*printSumFormula)(long power, FILE * out);

  /* To compute the sum for a specific power and the number of terms and obtain
   * CPU times
   * Parameters:
   *   power - desired power (IN)
   *   numTerms - number of terms (IN)
   *   result - sum computed (OUT)
   *   coeffInitTime - time in nanoseconds it took to initialize the
   *                   coefficients in the formula (OUT)
   *   sumTime - time to compute the sum (OUT)
   */
  void (*computeSumWithTimeStat)(long power, long numTerms, mpz_t result,
                                 long * coeffInitTime, long * sumTime);

  /* To compute the sum for a specific power and the number of terms
   * Parameters:
   *   power - desired power (IN)
   *   numTerms - number of terms (IN)
   *   result - sum computed (OUT)
   */
  void (*computeSum)(long power, long numTerms, mpz_t result);

  /* To compute the sum for a specific power and the number of terms using the
   * simple implementation (series summation)
   * Parameters:
   *   power - desired power (IN)
   *   numTerms - number of terms (IN)
   *   result - sum computed (OUT)
   */
  void (*computeSumUsingSeries)(long power, long numTerms, mpz_t result);

  /* To free the coefficients returned by getCoefficients()
   * Parameters:
   *   coeffs - coefficient array returned by getCoefficients() (IN)
   *   size - number of elements in the array (IN)
   */
  void (*freeCoefficients)(mpq_t * coeffs, long size);
} PowerSum;

void initPowerSum(PowerSum * psPtr);
void nCr(long n, long r, mpz_t result);
long computeCpuTime(struct timespec * before, struct timespec * after);

#endif
