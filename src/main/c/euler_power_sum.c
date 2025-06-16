#include <stdlib.h>

#include "power_sum.h"
#include "euler_power_sum.h"

static void printTerm(long start, long numTerms, FILE * out) {
  start++;
  for(long j = 0; j < numTerms; j++, start--) {
    if (start > 0) {
      fprintf(out, "(n + %ld)", start);
    } else if (start == 0) {
      fprintf(out, "n");
    } else {
      fprintf(out, "(n - %ld)", -start);
    }
  }
}

static void _freeEulerCoeffs(mpz_t * coeffs, long size) {
  long i;
  for (i = 0; i < size; i++) {
    mpz_clear(coeffs[i]);
  }
  free(coeffs);
}

/**
 * The Euler numbers of first kind are as defined below:
 * E(i, j) = 1 when j = 0
 * E(i, j) = (j + 1)E(n - 1, j) + (n - j)E(n - 1, j - 1)
 */
static void _getEulerCoefficients(long power, long maxN, mpz_t ** coeffList,
                                  long * numCoeffs) {

  if (power < 0) {
    *coeffList = NULL;
    *numCoeffs = 0;
    return;
  }

  mpz_t * coeffs;
  mpz_t temp;
  mpz_t E_i_1_j_1;
  long k;
  long j;

  long maxNumCoefficients = power + 1; /* Array is 0 based, so add 1 */

  mpz_init(temp);
  mpz_init(E_i_1_j_1);

  coeffs = (mpz_t *)malloc((maxNumCoefficients)*sizeof(mpz_t));
  mpz_init_set_ui(coeffs[0], 1UL);

  for (k = 1; k < maxNumCoefficients; k++) {
    mpz_init_set_ui(coeffs[k], 0UL);
  }

  for (long i = 1; i <= power; i++) {
    long oddPower = ((i & 1) == 1);
    long halfI = i >> 1;
    long halfLimit = oddPower ? halfI : (halfI - 1);
    long limit = halfLimit;
    /* No need to initialize more than maximum n since falling factorial in
     * other terms will be 0 */
    if (halfLimit > maxN) {
      limit = maxN;
    }

    for (j = 0; j <= limit; j++) {
      if (j == 0) {
        mpz_set(E_i_1_j_1, coeffs[j]);
        mpz_set_ui(coeffs[j], 1UL);
      } else {
        mpz_set(temp, coeffs[j]);
        mpz_mul_ui(temp, temp, (j + 1));
        mpz_mul_ui(E_i_1_j_1, E_i_1_j_1, (i - j));
        mpz_add(temp, temp, E_i_1_j_1);
        mpz_set(E_i_1_j_1, coeffs[j]);
        mpz_set(coeffs[j], temp);
      }
    }
    j = halfLimit;
    /* Since coefficients are mirror reflection w.r.t. the central point,
     * we can reverse copy whatever we initialized so far */
    if (oddPower) {
      for (k = 1; j + k < i; k++) {
        mpz_set(coeffs[j + k], coeffs[j - k]);
      }
    } else {
      for (k = 1; j + k < i; k++) {
        mpz_set(coeffs[j + k], coeffs[j - k + 1]);
      }
    }
    mpz_set_ui(coeffs[i], 0UL);
  }

  mpz_clear(temp);
  mpz_clear(E_i_1_j_1);
  *numCoeffs = maxNumCoefficients;
  *coeffList = coeffs;
}

void getEulerCoefficients(long power, mpq_t ** coeffList, long * numCoeffs) {
  mpz_t * coeffs;
  mpq_t * outCoeffs = NULL;
  _getEulerCoefficients(power, power, &coeffs, numCoeffs);
  if (*numCoeffs > 0) {
    outCoeffs = (mpq_t *)malloc(*numCoeffs*sizeof(mpq_t));
    for (long i = 0; i < *numCoeffs; i++) {
      mpq_init(outCoeffs[i]);
      mpq_set_z(outCoeffs[i], coeffs[i]);
    }
    _freeEulerCoeffs(coeffs, *numCoeffs);
  }
  *coeffList = outCoeffs;
}

void printEulerSumFormula(long power, FILE * out) {
  mpz_t * coeffs;
  long numCoeffs;
  mpz_t powerFactorial;
  long j;

  if (power < 0) {
    return;
  }

  mpz_init_set_ui(powerFactorial, 1UL);
  _getEulerCoefficients(power, power, &coeffs, &numCoeffs);
  /* Compute (power + 1)! */
  for (j = 2; j <= power + 1; j++) {
    mpz_mul_ui(powerFactorial, powerFactorial, j);
  }

  fprintf(out, "{ ");
  for (j = 0; j < numCoeffs; j++) {
    if (mpz_cmp_ui(coeffs[j], 1UL) > 0) {
      /* Coefficient is greater than 1 */
      gmp_fprintf(out, "%Zd", coeffs[j]);
    }
    if (mpz_cmp_ui(coeffs[j], 0UL) > 0) {
      printTerm(j, power + 1, out);
    }
    if (j < (numCoeffs - 2)) {
      fprintf(out, " + ");
    }
  }
  fprintf(out, " }");
  if (mpz_cmp_ui(powerFactorial, 1UL) != 0) {
    gmp_fprintf(out, "/%Zd", powerFactorial);
  }
  fprintf(out, "\n");
  _freeEulerCoeffs(coeffs, numCoeffs);
  mpz_clear(powerFactorial);
}

void computeEulerSumWithTimeStat(long power, long n, mpz_t result,
                                 long * coeffInitTime, long * sumTime) {
  struct timespec before;
  struct timespec after;
  mpz_t * coeffs;
  long numCoeffs;

  mpz_set_ui(result, 0UL);
  if (power < 0 || n < 0) {
    *coeffInitTime = 0;
    *sumTime = 0;
    return;
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  _getEulerCoefficients(power, n, &coeffs, &numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *coeffInitTime = computeCpuTime(&before, &after);

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  long j;
  mpz_t temp;
  mpz_t temp1;
  unsigned long firstTime = 1;

  mpz_init(temp);
  mpz_init(temp1);
  for (j = 0; j < numCoeffs; j++) {
    if (n + j >= power) {
      if (firstTime) {
        firstTime = 0;
        nCr(n + j + 1, power + 1, temp);
      } else {
        mpz_mul_ui(temp, temp, n + j + 1);
        mpz_divexact_ui(temp, temp, n + j - power);
      }
      mpz_mul(temp1, coeffs[j], temp);
      mpz_add(result, result, temp1);
    }
  }
  mpz_clear(temp1);
  mpz_clear(temp);
  _freeEulerCoeffs(coeffs, numCoeffs);

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *sumTime = computeCpuTime(&before, &after);
}

void computeEulerSum(long power, long n, mpz_t result) {
  long cTime;
  long sTime;

  computeEulerSumWithTimeStat(power, n, result, &cTime, &sTime);
}

void freeEulerCoeffs(mpq_t * coeffs, long size) {
  long i;
  for (i = 0; i < size; i++) {
    mpq_clear(coeffs[i]);
  }
  free(coeffs);
}

void initEulerPowerSum(PowerSum * psPtr) {
  initPowerSum(psPtr);
  psPtr->getCoefficients = getEulerCoefficients;
  psPtr->printSumFormula = printEulerSumFormula;
  psPtr->computeSumWithTimeStat = computeEulerSumWithTimeStat;
  psPtr->computeSum = computeEulerSum;
  psPtr->freeCoefficients = freeEulerCoeffs;
}
