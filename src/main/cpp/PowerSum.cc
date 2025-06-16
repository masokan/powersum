#include "PowerSum.h"

mpz_class PowerSum :: computeSumUsingSeries(long power, long n) {
  mpz_class sum = 0;
  if (power < 0 || n < 0) {
    return sum;
  }

  for (long count = 0; count <= n; count++) {
    mpz_class base = count;
    mpz_class result;
    mpz_pow_ui(result.get_mpz_t(), base.get_mpz_t(), (unsigned long)power);
    sum += result;
  }
  return sum;
}

mpz_class PowerSum :: computeSum(long power, long n) {
  vector<long> stat;
  return(computeSumWithTimeStat(power, n, stat));
}

long PowerSum :: computeCpuTime(struct timespec & before,
                                struct timespec & after) {
  return((after.tv_sec - before.tv_sec)*1000000000L
         + after.tv_nsec - before.tv_nsec);
}

mpz_class PowerSum :: nCr(long n, long r) {
  long num = n;
  long i;

  if (r > n) {
    return 0;
  }
  if (r > n/2) {
    r = n - r;
  }
  if (r == 0) {
    return 1;
  }
  mpz_class result = num;
  for (i = 2; i <= r; i++) {
    num--;
    result *= num;
    result /= i;
  }
  return(result);
}
