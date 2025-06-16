#include <stdlib.h>

#include "power_sum.h"
#include "bernoulli_power_sum.h"

static void computeNextCoefficient(mpq_t * current, long m, mpq_t nextCoeff) {
  mpz_t binom;
  mpq_t temp;

  mpz_init_set_ui(binom, 1UL);
  mpq_set_ui(nextCoeff, 0UL, 1UL);
  mpq_init(temp);

  for (long k = 0; k < m; k++) {
    if (((k & 1) == 0) || k == 1) {
      /* Even coefficient or the first odd one */
      mpq_set_z(temp, binom);
      mpq_mul(temp, temp, current[k]);
      mpq_add(nextCoeff, nextCoeff, temp);
    }
    mpz_mul_ui(binom, binom, (unsigned long)(m + 1 - k));
    mpz_divexact_ui(binom, binom, (unsigned long)(k + 1));
  }
  mpq_neg(nextCoeff, nextCoeff);
  mpq_set_z(temp, binom);
  mpq_div(nextCoeff, nextCoeff, temp);

  mpq_clear(temp);
  mpz_clear(binom);
}

void freeBernoulliCoeffs(mpq_t * coeffs, long size) {
  long i;
  for (i = 0; i < size; i++) {
    mpq_clear(coeffs[i]);
  }
  free(coeffs);
}

/* Bernoulli numbers are defined by the following recurrence relation:
 * B(0) = 1
 * B(m) = -(Binom((m + 1), 0)B(0) + Binom((m + 1), 1)B(1)
 *          + ... + Binom((m + 1), (m - 1))B(m - 1))
 * where Binom(i, j) is the binomial coefficient which evaluates to:
 * i!/{(i - j)!j!}
 */
void getBernoulliCoefficients(long power, mpq_t ** coeffList,
                              long * numCoeffs) {
  if (power < 0) {
    *coeffList = NULL;
    *numCoeffs = 0;
    return;
  }

  mpq_t * coeffs = (mpq_t *)malloc((power + 1)*sizeof(mpq_t));
  long k;

  for (k = 0; k <= power; k++) {
    mpq_init(coeffs[k]);
  }

  /* Initialize B(0) */
  mpq_set_ui(coeffs[0], 1UL, 1UL);
  if (power > 0) {
    /* Initialize B(1) */
    mpq_set_si(coeffs[1], -1, 2UL);

    for (long i = 2; i <= power; i++) {
      if (i & 1) {
        /* Odd coefficients above 1 are 0s */
        mpq_set_ui(coeffs[i], 0UL, 1UL);
      } else {
        computeNextCoefficient(coeffs, i, coeffs[i]);
      }
    }
  }
  *numCoeffs = power + 1;
  *coeffList = coeffs;
}

void printBernoulliSumFormula(long power, FILE * out) {

  if (power < 0) {
    return;
  }

  mpq_t * coeffs;
  long numCoeffs;
  mpz_t binom;
  long binomN = power + 1;
  long binomR = 1;
  long p = power + 1;
  mpz_t one;

  mpz_init_set_ui(one, 1);
  mpz_init_set_ui(binom, 1);
  getBernoulliCoefficients(power, &coeffs, &numCoeffs);
  fprintf(out, "{ ");
  for (long i = 0; i < numCoeffs; i++) {
    if ((i & 1) == 0 || i == 1) {
      if (i != 0) {
        fprintf(out, " + ");
      }
      if (mpq_cmp_z(coeffs[i], one) != 0) {
        gmp_fprintf(out, "(%Qd)", coeffs[i]);
      }
      if (mpz_cmp(binom, one) != 0) {
        gmp_fprintf(out, "%Zd", binom);
      }
      fprintf(out, "(n + 1)");
      if (p != 1) {
        fprintf(out, "^%ld", p);
      }
    }
    p--;
    mpz_mul_ui(binom, binom, binomN);
    mpz_divexact_ui(binom, binom, binomR);
    binomN--;
    binomR++;
  }
  fprintf(out, " }");
  if (power > 0) {
    fprintf(out, "/%ld", (power + 1));
  }
  fprintf(out, "\n");
  freeBernoulliCoeffs(coeffs, numCoeffs);
  mpz_clear(one);
  mpz_clear(binom);
}

void computeBernoulliSumWithTimeStat(long power, long n, mpz_t result,
                                     long * coeffInitTime, long * sumTime) {
  mpq_t * coeffs;
  mpq_t sum;
  mpq_t temp;
  mpz_t binom;
  mpz_t pow;
  unsigned long binomN = power + 1;
  unsigned long binomR = 1;
  long numCoeffs;
  mpq_t binomq;
  mpq_t powq;

  struct timespec before;
  struct timespec after;

  if (power < 0 || n < 0) {
    mpz_set_ui(result, 0UL);
    *coeffInitTime = 0;
    *sumTime = 0;
    return;
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  getBernoulliCoefficients(power, &coeffs, &numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *coeffInitTime = computeCpuTime(&before, &after);

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  mpq_init(sum);
  mpq_set_ui(sum, 0UL, 1UL);
  mpq_init(temp);
  mpz_init_set_ui(binom, binomN);
  long nPlusOne = n + 1;
  mpz_init_set_ui(pow, nPlusOne);
  mpq_init(binomq);
  mpq_init(powq);

  /* We compute the terms in reverse order to avoid an additional division
     operation in the loop */
  for (long i = numCoeffs - 1; i >= 0 ; i--) {
    if ((i & 1) == 0 || i == 1) {
      mpq_set(temp, coeffs[i]);
      mpq_set_z(binomq, binom);
      mpq_mul(temp, temp, binomq);
      mpq_set_z(powq, pow);
      mpq_mul(temp, temp, powq);
      mpq_add(sum, sum, temp);
    }
    mpz_mul_ui(pow, pow, nPlusOne);
    binomN--;
    binomR++;
    mpz_mul_ui(binom, binom, binomN);
    mpz_divexact_ui(binom, binom, binomR);
  }
  mpq_set_ui(powq, (unsigned long)(power + 1), 1UL);
  mpq_div(sum, sum, powq);
  mpq_canonicalize(sum);
  freeBernoulliCoeffs(coeffs, numCoeffs);
  mpz_clear(pow);
  mpz_clear(binom);
  mpq_clear(temp);
  mpq_clear(powq);
  mpq_clear(binomq);
  mpq_get_num(result, sum);
  mpq_clear(sum);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *sumTime = computeCpuTime(&before, &after);
}

void computeBernoulliSum(long power, long n, mpz_t result) {
  long cTime;
  long sTime;

  computeBernoulliSumWithTimeStat(power, n, result, &cTime, &sTime);
}

void initBernoulliPowerSum(PowerSum * psPtr) {
  initPowerSum(psPtr);
  psPtr->getCoefficients = getBernoulliCoefficients;
  psPtr->printSumFormula = printBernoulliSumFormula;
  psPtr->computeSumWithTimeStat = computeBernoulliSumWithTimeStat;
  psPtr->computeSum = computeBernoulliSum;
  psPtr->freeCoefficients = freeBernoulliCoeffs;
}
