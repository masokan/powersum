#include <iostream>

using std::endl;

#include "EulerPowerSum.h"

EulerPowerSum :: EulerPowerSum()
                  : PowerSum() {
}

vector<mpq_class> EulerPowerSum :: getCoefficients(long power) {
  vector<mpq_class> outCoeffs;
  vector<mpz_class> coeffs = getCoefficients(power, power);

  for (size_t i = 0; i < coeffs.size(); i++) {
    outCoeffs.push_back(mpq_class(coeffs[i]));
  }
  return outCoeffs;
}

void EulerPowerSum :: printSumFormula(long power, ostream &out) {
  if (power < 0) {
    return;
  }

  vector<mpz_class> coeffs = getCoefficients(power, power);
  mpz_class powerFactorial;

  // Compute (power + 1)!
  powerFactorial = 1;
  for (long i = 2; i <= power + 1; i++) {
    powerFactorial *= i;
  }

  out << "{ ";
  long numCoeffs = (long)coeffs.size();
  for (long j = 0; j < numCoeffs; j++) {
    if (coeffs[j] > 1) {
      out << coeffs[j];
    }
    if (coeffs[j] > 0) {
      printTerm(j, power + 1, out);
    }
    if (j < (numCoeffs - 2)) {
      out << " + ";
    }
  }
  out << " }";
  if (powerFactorial != 1) {
    out << '/' << powerFactorial;
  }
  out << endl;
}

mpz_class EulerPowerSum :: computeSumWithTimeStat(long power, long n, 
                                                  vector<long> & stat) {
  struct timespec before;
  struct timespec after;

  mpz_class sum = 0;
  stat.clear();

  if (power < 0 || n < 0) {
    stat.push_back(0);
    stat.push_back(0);
    return sum;
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  vector<mpz_class> coeffs = getCoefficients(power, n);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  mpz_class temp;

  bool firstTime = true;
  long numCoeffs = (long)coeffs.size();
  for (long j = 0; j < numCoeffs; j++) {
    if (n + j >= power) {
      if (firstTime) {
        firstTime = false;
        temp = nCr(n + j + 1, power + 1);
      } else {
        temp *= (n + j + 1);
        temp /= (n + j - power);
      }
      sum += coeffs[j]*temp;
    }
  }

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));
  return sum;
}

EulerPowerSum :: ~EulerPowerSum() {
}

/**
 * The Euler numbers of first kind are as defined below:
 * E(i, j) = 1 when j = 0
 * E(i, j) = (j + 1)E(n - 1, j) + (n - j)E(n - 1, j - 1)
 */
vector<mpz_class> EulerPowerSum :: getCoefficients(long power,
                                                   long maxNumCoefficients) {
  vector<mpz_class> coeffs;

  if (power < 0) {
    return coeffs;
  }

  mpz_class temp;
  mpz_class E_i_1_j_1;

  long k;
  long j;

  coeffs.push_back(1);
  for (k = 1; k < power + 1; k++) {
    coeffs.push_back(0);
  }

  for (long i = 1; i <= power; i++) {
    bool oddPower = ((i & 1) == 1);
    long halfI = i >> 1;
    long halfLimit = oddPower ? halfI : (halfI - 1);
    long limit = halfLimit;
    // No need to initialize more than maximum n since falling factorial in
    // other terms will be 0
    if (halfLimit > maxNumCoefficients) {
      limit = maxNumCoefficients;
    }

    for (j = 0; j <= limit; j++) {
      if (j == 0) {
        E_i_1_j_1 = coeffs[j];
        coeffs[j] = 1;
      } else {
        temp = (j + 1)*coeffs[j] + (i - j)*E_i_1_j_1;
        E_i_1_j_1 = coeffs[j];
        coeffs[j] = temp;
      }
    }
    j = halfLimit;
    // Since coefficients are mirror reflection w.r.t. the central point,
    // we can reverse copy whatever we initialized so far
    if (oddPower) {
      for (k = 1; j + k < i; k++) {
        coeffs[j + k] = coeffs[j - k];
      }
    } else {
      for (k = 1; j + k < i; k++) {
        coeffs[j + k] = coeffs[j - k + 1];
      }
    }
    coeffs[i] = 0;
  }
  return coeffs;
}

void EulerPowerSum :: printTerm(long start, long numTerms, ostream & out) {
  start++;
  for(long j = 0; j < numTerms; j++, start--) {
    if (start > 0) {
      out << "(n + " << start << ')';
    } else if (start == 0) {
      out << 'n';
    } else {
      out << "(n - " << -start << ')';
    }
  }
}
