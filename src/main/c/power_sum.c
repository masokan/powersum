#include <stdlib.h>

#include "power_sum.h"

void computeSumUsingSeries(long power, long n, mpz_t result) {
  long count;
  mpz_t oneTerm;
  mpz_set_ui(result, 0UL);

  if (power < 0 || n < 0) {
    return;
  }

  mpz_init(oneTerm);
  for (count = 0; count <= n; count++) {
    mpz_ui_pow_ui(oneTerm, (unsigned long)count, (unsigned long)power);
    mpz_add(result, result, oneTerm);
  }
  mpz_clear(oneTerm);
}

void initPowerSum(PowerSum * psPtr) {
  psPtr->getCoefficients = NULL;
  psPtr->printSumFormula = NULL;
  psPtr->computeSumWithTimeStat = NULL;
  psPtr->computeSum = NULL;
  psPtr->computeSumUsingSeries = computeSumUsingSeries;
  psPtr->freeCoefficients = NULL;
}

void nCr(long n, long r, mpz_t result) {
  long num = n;
  long i;

  if (r > n) {
    mpz_set_ui(result, 0UL);
    return;
  }
  if (r > n/2) {
    r = n - r;
  }
  if (r == 0) {
    mpz_set_ui(result, 1UL);
    return;
  }
  mpz_set_ui(result, num);

  for (i = 2; i <= r; i++) {
    num--;
    mpz_mul_ui(result, result, num);
    mpz_divexact_ui(result, result, i);
  }
}

long computeCpuTime(struct timespec * before, struct timespec * after) {
  return((after->tv_sec - before->tv_sec)*1000000000L
         + after->tv_nsec - before->tv_nsec);
}
