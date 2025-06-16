#include <iostream>

using std::endl;

#include "BernoulliPowerSum.h"

BernoulliPowerSum :: BernoulliPowerSum()
                   : PowerSum() {
}

/* Bernoulli numbers are defined by the following recurrence relation:
 * B(0) = 1
 * B(m) = -(Binom((m + 1), 0)B(0) + Binom((m + 1), 1)B(1)
 *          + ... + Binom((m + 1), (m - 1))B(m - 1))
 * where Binom(i, j) is the binomial coefficient which evaluates to:
 * i!/{(i - j)!j!}
 */
vector<mpq_class> BernoulliPowerSum :: getCoefficients(long power) {
  vector<mpq_class> coeffs;

  if (power < 0) {
    return coeffs;
  }

  // Initialize B(0)
  coeffs.push_back(1);
  if (power > 0) {
    // Initialize B(1)
    coeffs.push_back(mpq_class(-1, 2));

    for (long i = 2; i <= power; i++) {
      if (i & 1) {
        // Odd coefficients above 1 are 0s
        coeffs.push_back(0);
      } else {
        coeffs.push_back(computeNextCoefficient(coeffs, i));
      }
    }
  }
  return coeffs;
}

void BernoulliPowerSum :: printSumFormula(long power, ostream &out) {
  if (power < 0) {
    return;
  }

  vector<mpq_class> coeffs = getCoefficients(power);
  mpz_class binom = 1;
  long binomN = power + 1;
  long binomR = 1;
  long p = power + 1;

  out << "{ ";
  for (size_t i = 0; i < coeffs.size(); i++) {
    if ((i & 1) == 0 || i == 1) {
      if (i != 0) {
        out << " + ";
      }
      if (coeffs[i] != 1) {
        out << '(' << coeffs[i] << ')';
      }
      if (binom != 1) {
        out << binom;
      }
      out << "(n + 1)";
      if (p != 1) {
        out << '^' << p;
      }
    }
    p--;
    binom = binom*binomN/binomR;
    binomN--;
    binomR++;
  }
  out << " }";
  if (power > 0) {
    out << '/' << (power + 1);
  }
  out << endl;
}

mpz_class BernoulliPowerSum :: computeSumWithTimeStat(long power, long n, 
                                                      vector<long> & stat) {
  struct timespec before;
  struct timespec after;

  mpq_class sum = 0;
  stat.clear();

  if (power < 0 || n < 0) {
    stat.push_back(0);
    stat.push_back(0);
    return sum.get_num();
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  vector<mpq_class> coeffs = getCoefficients(power);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  mpz_class pow = n + 1;
  long binomN = power + 1;
  long binomR = 1;
  mpz_class binom = binomN;

  // We compute the terms in reverse order to avoid an additional division
  // operation in the loop
  for (long i = coeffs.size() - 1; i >= 0; i--) {
    if ((i & 1) == 0 || i == 1) {
      sum += binom*coeffs[i]*pow;
    }
    pow *= (n + 1);
    binomN--;
    binomR++;
    binom = binom*binomN/binomR;
  }
  sum /= (power + 1);
  sum.canonicalize();
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));
  return sum.get_num();
}

BernoulliPowerSum :: ~BernoulliPowerSum() {
}

mpq_class BernoulliPowerSum :: computeNextCoefficient(
                                       vector<mpq_class> & current, long m) {
  mpz_class binom = 1;
  mpq_class coeff = 0;

  for (long k = 0; k < m; k++) {
    if (((k & 1) == 0) || k == 1) {
      // Even coefficient or the first odd one
      coeff += current[k]*binom;
    }
    binom *= (m + 1 - k);
    binom /= (k + 1);
  }
  coeff = -coeff/binom;
  coeff.canonicalize();
  return coeff;
}
