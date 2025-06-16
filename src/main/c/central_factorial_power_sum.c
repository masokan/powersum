#include <stdlib.h>

#include "power_sum.h"
#include "central_factorial_power_sum.h"

static void printFallingFactorial(long start, long numTerms, FILE * out) {
  for(long i = 0; i < numTerms; i++) {
    if (i == 0) {
      fprintf(out, "(n + %ld)", start);
    } else {
      long diff = start - i;
      if (diff == 0) {
        fprintf(out, "n");
      } else if (diff > 0) {
        fprintf(out, "(n + %ld)", diff);
      } else {
        fprintf(out, "(n - %ld)", -diff);
      }
    }
  }
}

static void _freeCentralFactorialCoeffs(mpz_t * coeffs, long size) {
  for (long i = 0; i < size; i++) {
    mpz_clear(coeffs[i]);
  }
  free(coeffs);
}

/**
 * The Central Factorial Numbers of second kind are as defined below:
 * T(2m, 2k) = k*k*T(2m - 2, 2k) + T(2m - 2, 2k - 2)
 * T(2m, 2m) = 1
 * We need only even coefficients.
 */
static void _getCentralFactorialCoefficients(long power, long maxN,
                                             mpz_t ** coeffList,
                                             long * numCoeffs) {
  mpz_t * coeffs;
  mpz_t T_2_2;
  mpz_t temp;
  long i;
  long k;

  if (power < 0) {
    *coeffList = NULL;
    *numCoeffs = 0;
    return;
  }

  /* Get half of power.  For odd power, round up half power */
  long m = (power >> 1) + (power & 1);

  long maxNumCoefficients = m + 1; /* Array is 0 based, so add 1 */

  /* We do not compute coefficients more than what is needed.  Example: power
   * is 1000 and the number of terms in the series is only 3 like
   * (1^1000 + 2^1000 + 3^1000).  In this case, falling factorials in many
   * terms will be 0.  There is no point in computing the coefficients for
   * these terms.
   */
  if (maxN < maxNumCoefficients) {
    maxNumCoefficients = maxN + 1;
  }

  mpz_init(temp);
  mpz_init(T_2_2);

  coeffs = (mpz_t *)malloc((maxNumCoefficients)*sizeof(mpz_t));
  for (k = 0; k < maxNumCoefficients; k++) {
    mpz_init(coeffs[k]);
  }

  for (i = 0; i <= m; i++) {
    for (k = 0; k < maxNumCoefficients; k++) {
      if (i == k) {
        mpz_set_ui(coeffs[k], 1UL);
      } else if (i > 0 && k > 0) {
        /* Calculate T(2*i, 2*k) = T(2*i-2,2*k-2) + k*k*T(2*i-2, 2*k)
         * Use the previous values in coeffs[k] for all the computation before
         * setting the new value.
         */
        mpz_mul_ui(temp, coeffs[k], k);
        mpz_mul_ui(temp, temp, k);
        mpz_add(temp, temp, T_2_2);
        mpz_set(T_2_2, coeffs[k]);
        mpz_set(coeffs[k], temp);
      } else {
        mpz_set_ui(coeffs[k], 0UL);
      }
    }
    mpz_set(T_2_2, coeffs[0]);
  }

  mpz_clear(T_2_2);
  mpz_clear(temp);
  *numCoeffs = maxNumCoefficients;
  *coeffList = coeffs;
}

void getCentralFactorialCoefficients(long power, mpq_t ** coeffList,
                                     long * numCoeffs) {
  mpz_t * coeffs;
  mpq_t * outCoeffs = NULL;
  _getCentralFactorialCoefficients(power, power, &coeffs, numCoeffs);
  if (*numCoeffs > 0) {
    outCoeffs = (mpq_t *)malloc(*numCoeffs*sizeof(mpq_t));
    for (long i = 0; i < *numCoeffs; i++) {
      mpq_init(outCoeffs[i]);
      mpq_set_z(outCoeffs[i], coeffs[i]);
    }
    _freeCentralFactorialCoeffs(coeffs, *numCoeffs);
  }
  *coeffList = outCoeffs;
}

void printCentralFactorialSumFormula(long power, FILE * out) {
  mpz_t * coeffs;
  long numCoeffs;

  if (power < 0) {
    return;
  }

  long evenPower = ((power & 1) == 0);
  _getCentralFactorialCoefficients(power, power, &coeffs, &numCoeffs);
  if (power > 0) {
    /* Coefficient 0 is always 0 - so we start the loop at index 1 */
    for (long k = 1; k < numCoeffs; k++) {
      if (mpz_cmp_ui(coeffs[k], 1UL) > 0) {
        /* Coefficient is greater than 1 */
        gmp_fprintf(out, "%Zd", coeffs[k]);
      }
      if (evenPower) {
        fprintf(out, "(2n + 1)");
      }
      printFallingFactorial(k, 2*k, out);
      if (evenPower) {
        fprintf(out, "/%ld", 2*(2*k + 1));
      } else {
        fprintf(out, "/%ld", 2*k);
      }
      if (k != numCoeffs - 1) {
        fprintf(out, " + ");
      }
    }
  } else {
    /* Special case for 0 */
    fprintf(out, "(n + 1)");
  }
  fprintf(out, "\n");
  _freeCentralFactorialCoeffs(coeffs, numCoeffs);
}

void computeCentralFactorialSumWithTimeStat(long power, long n, mpz_t result,
                                        long * coeffInitTime, long * sumTime) {
  struct timespec before;
  struct timespec after;
  mpz_t * coeffs;
  long numCoeffs;

  if (power < 0 || n < 0) {
    *coeffInitTime = 0;
    *sumTime = 0;
    return;
  }

  long evenPower = ((power & 1) == 0);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  _getCentralFactorialCoefficients(power, n, &coeffs, &numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *coeffInitTime = computeCpuTime(&before, &after);

  if (power > 0) {
    mpz_set_ui(result, 0UL);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    mpz_t fallingFactorial;
    mpz_t temp;

    mpz_init_set_ui(fallingFactorial, 1);
    mpz_init(temp);

    /* Coefficient 0 is always 0 - so we start the loop at index 1 */
    for (long k = 1; k < numCoeffs; k++) {
      mpz_mul_ui(fallingFactorial, fallingFactorial, n + k);
      mpz_mul_ui(fallingFactorial, fallingFactorial, n - k +1);
      mpz_mul(temp, coeffs[k], fallingFactorial);
      if (evenPower) {
        mpz_mul_ui(temp, temp, (2*n + 1));
        mpz_divexact_ui(temp, temp, (2*(2*k + 1)));
      } else {
        mpz_divexact_ui(temp, temp, 2*k);
      }
      mpz_add(result, result, temp);
    }
    mpz_clear(fallingFactorial);
    mpz_clear(temp);
  } else {
    /* Special case - not handled by the formula */
    mpz_set_ui(result, n + 1);
  }
  _freeCentralFactorialCoeffs(coeffs, numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *sumTime = computeCpuTime(&before, &after);
}

void computeCentralFactorialSum(long power, long n, mpz_t result) {
  long cTime;
  long sTime;

  computeCentralFactorialSumWithTimeStat(power, n, result, &cTime, &sTime);
}

void freeCentralFactorialCoeffs(mpq_t * coeffs, long size) {
  for (long i = 0; i < size; i++) {
    mpq_clear(coeffs[i]);
  }
  free(coeffs);
}

void initCentralFactorialPowerSum(PowerSum * psPtr) {
  initPowerSum(psPtr);
  psPtr->getCoefficients = getCentralFactorialCoefficients;
  psPtr->printSumFormula = printCentralFactorialSumFormula;
  psPtr->computeSumWithTimeStat = computeCentralFactorialSumWithTimeStat;
  psPtr->computeSum = computeCentralFactorialSum;
  psPtr->freeCoefficients = freeCentralFactorialCoeffs;
}
