#include <stdlib.h>

#include "stirling_power_sum.h"

/**
  * The coefficients are Stirling numbers of second kind which are defined as
  * below:
  * S(0, 0) = 1
  * S(m, 0) = 0 for all m > 0
  * S(m, j) = 0 for all j > m
  * S(m, j) = S(m-1, j-1) + j*S(m-1, j)
  * m is the power of the series (m >= 0)
  * j is the index of the term in the sum formula (j >= 0)
  */
static void _getStirlingCoefficients(long power, long maxNumCoefficients,
                                     mpz_t ** coeffList, long * numCoeffs) {
  if (power < 0) {
    *coeffList = NULL;
    *numCoeffs = 0;
    return;
  }

  mpz_t * coeffs = (mpz_t *)malloc((power + 1)*sizeof(mpz_t));
  *numCoeffs = power + 1;
  mpz_t Sm_1j_1;
  mpz_t Smj;
  long currentPower;
  long term;

  mpz_init_set_ui(Sm_1j_1, 1UL); /* S(0, 0) */
  mpz_init_set_ui(Smj, 0UL);
  if (power > 0) {
    maxNumCoefficients--;
    /* Initialize all coefficients to 0 */
    for (currentPower = 0;
         currentPower <= power;
         currentPower++) {
      mpz_init_set_ui(coeffs[currentPower], 0UL);
    }

    for (currentPower = 1; currentPower <= power; currentPower++) {
      long termMax
       = (currentPower > maxNumCoefficients)? maxNumCoefficients : currentPower;
      /* Initialize coefficients for the current power */
      for(term = 1; term <= termMax; term++) {
        mpz_set_ui(Smj, term);
        mpz_mul(Smj, Smj, coeffs[term]);
        mpz_add(Smj, Smj, Sm_1j_1);
        mpz_set(Sm_1j_1, coeffs[term]);
        mpz_set(coeffs[term], Smj);
      }
      mpz_set(Sm_1j_1, coeffs[0]);
    }
  } else {
    mpz_init_set_ui(coeffs[0], 1UL);
  }
  mpz_clear(Sm_1j_1);
  mpz_clear(Smj);
  *coeffList = coeffs;
}

static void printStirlingFactors(long term, FILE * out) {
  long i;
  for(i = 0; i <= term; i++) {
    switch(i) {
      case 0:
        fprintf(out, "(n + 1)");
        break;
      case 1:
        fprintf(out, "n");
        break;
      default:
        fprintf(out, "(n - %ld)", (i - 1));
        break;
    }
  }
}

static void _freeStirlingCoeffs(mpz_t * coeffs, long size) {
  for (long i = 0; i < size; i++) {
    mpz_clear(coeffs[i]);
  }
  free(coeffs);
}

void getStirlingCoefficients(long power, mpq_t ** coeffList, long * numCoeffs) {
  mpz_t * coeffs;
  mpq_t * outCoeffs = NULL;
  _getStirlingCoefficients(power, power + 1, &coeffs, numCoeffs);
  if (*numCoeffs > 0) {
    outCoeffs = (mpq_t *)malloc(*numCoeffs*sizeof(mpq_t));
    for (long i = 0; i < *numCoeffs; i++) {
      mpq_init(outCoeffs[i]);
      mpq_set_z(outCoeffs[i], coeffs[i]);
    }
    _freeStirlingCoeffs(coeffs, *numCoeffs);
  }
  *coeffList = outCoeffs;
}

void freeStirlingCoeffs(mpq_t * coeffs, long size) {
  for (long i = 0; i < size; i++) {
    mpq_clear(coeffs[i]);
  }
  free(coeffs);
}

void printStirlingSumFormula(long power, FILE * out) {
  long t;
  mpz_t * coeffs;
  long numCoeffs;

  if (power < 0) {
    return;
  }

  _getStirlingCoefficients(power, power + 1, &coeffs, &numCoeffs);
  unsigned int firstTime = 1;
  for (t = 0; t <= numCoeffs - 1; t++) {
    if (mpz_cmp_ui(coeffs[t], 0UL) != 0) {
      if (firstTime == 1) {
        fprintf(out, "   ");
        firstTime = 0;
      } else {
        fprintf(out, " + ");
      }
      if (mpz_cmp_ui(coeffs[t], 1UL) > 0) {
        /* Coefficient is greater than 1 */
        gmp_fprintf(out, "%Zd", coeffs[t]);
      }
      printStirlingFactors(t, out);
      if (t > 0) {
        fprintf(out, "/%ld", (t + 1));
      }
    }
  }
  fprintf(out, "\n");
  _freeStirlingCoeffs(coeffs, numCoeffs);
}

void computeStirlingSumWithTimeStat(long power, long n, mpz_t result,
                                    long * coeffInitTime, long * sumTime) {
  struct timespec before;
  struct timespec after;
  mpz_t fallingFactorial;
  mpz_t factor;
  mpz_t * coeffs;
  long numCoeffs;
  mpz_set_ui(result, 0UL);

  if (power < 0 || n < 0) {
    *coeffInitTime = 0;
    *sumTime = 0;
    return;
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  _getStirlingCoefficients(power, n + 1, &coeffs, &numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *coeffInitTime = computeCpuTime(&before, &after);

  long numTermsToCompute;
  long t;

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  numTermsToCompute = (power > n)? n : power;
  numTermsToCompute++;
  mpz_init(factor);
  mpz_init_set_ui(fallingFactorial, n + 1);
  for (t = 0; t < numTermsToCompute; t++) {
    /* We know that the following division is exact and we do this first before
     * multiplication to avoid generating a large intermediate value
     */
    mpz_divexact_ui(factor, fallingFactorial, t + 1);
    mpz_mul(factor, factor, coeffs[t]);
    mpz_add(result, result, factor);
    mpz_mul_ui(fallingFactorial, fallingFactorial, n - t);
  }
  mpz_clear(fallingFactorial);
  mpz_clear(factor);
  _freeStirlingCoeffs(coeffs, numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *sumTime = computeCpuTime(&before, &after);
}

void computeStirlingSum(long power, long n, mpz_t result) {
  long cTime;
  long sTime;

  computeStirlingSumWithTimeStat(power, n, result, &cTime, &sTime);
}

void initStirlingPowerSum(PowerSum * psPtr) {
  initPowerSum(psPtr);
  psPtr->getCoefficients = getStirlingCoefficients;
  psPtr->printSumFormula = printStirlingSumFormula;
  psPtr->computeSumWithTimeStat = computeStirlingSumWithTimeStat;
  psPtr->computeSum = computeStirlingSum;
  psPtr->freeCoefficients = freeStirlingCoeffs;
}
