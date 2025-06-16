#include <iostream>

using std::endl;

#include "CentralFactorialPowerSum.h"

CentralFactorialPowerSum :: CentralFactorialPowerSum()
                          : PowerSum() {
}

vector<mpq_class> CentralFactorialPowerSum :: getCoefficients(long power) {
  vector<mpq_class> outCoeffs;
  vector<mpz_class> coeffs = getCoefficients(power, power);

  for (size_t i = 0; i < coeffs.size(); i++) {
    outCoeffs.push_back(mpq_class(coeffs[i]));
  }
  return outCoeffs;
}

void CentralFactorialPowerSum :: printSumFormula(long power, ostream &out) {
  if (power < 0) {
    return;
  }

  bool evenPower = ((power & 1) == 0);

  if (power > 0) {
    vector<mpz_class>  coeffs = getCoefficients(power, power);
    long numCoeffs = coeffs.size();
    // Coefficient 0 is always 0 - so we start the loop at index 1
    for (long k = 1; k < numCoeffs; k++) {
      mpz_class coeff = coeffs[k];
      if (coeff != 1) {
        out << coeff;
      }
      if (evenPower) {
        out << "(2n + 1)";
      }
      printFallingFactorial(k, 2*k, out);
      if (evenPower) {
        out << '/' << 2*(2*k + 1);
      } else {
        out << '/' << 2*k;
      }
      if (k < numCoeffs - 1) {
        out << " + ";
      }
    }
  } else {
    // Special case for 0
    out << "(n + 1)";
  }
  out << endl;
}

mpz_class CentralFactorialPowerSum :: computeSumWithTimeStat(long power,
                                             long n, vector<long> & stat) {
  struct timespec before;
  struct timespec after;

  mpz_class sum = 0;
  stat.clear();

  if (power < 0 || n < 0) {
    stat.push_back(0);
    stat.push_back(0);
    return sum;
  }

  if (power > 0) {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    vector<mpz_class> coeffs = getCoefficients(power, n);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
    stat.push_back(computeCpuTime(before, after));

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    bool evenPower = ((power & 1) == 0);
    long numCoeffs = coeffs.size();
    mpz_class fallingFactorial = 1;

    // Coefficient 0 is always 0 - so we start the loop at index 1
    for (long k = 1; k < numCoeffs; k++) {
      fallingFactorial *= (n + k)*(n - k + 1);
      if (evenPower) {
        sum += (coeffs[k]*fallingFactorial*(2*n + 1))/(2*(2*k + 1));
      } else {
        sum += coeffs[k]*fallingFactorial/(2*k);
      }
    }
  } else {
    // Special case - not handled by the formula
    stat.push_back(0);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    sum = n + 1;
  }
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));
  return sum;

}

CentralFactorialPowerSum :: ~CentralFactorialPowerSum() {
}

/**
 * The Central Factorial Numbers of second kind are as defined below:
 * T(2m, 2k) = k*k*T(2m - 2, 2k) + T(2m - 2, 2k - 2)
 * T(2m, 2m) = 1
 * We need only even coefficients.
 */
vector<mpz_class> CentralFactorialPowerSum :: getCoefficients(long power,
                                                              long maxN) {
  vector<mpz_class> coeffs;

  if (power < 0) {
    return coeffs;
  }

  mpz_class T_2_2 = 0;
  mpz_class temp;
  long i;
  long k;

  // Get half of power.  For odd power, round up half power
  long m = (power >> 1) + (power & 1);

  long maxNumCoefficients = m + 1; // Array is 0 based, so add 1

  // We do not compute coefficients more than what is needed.  Example: power
  // is 1000 and the number of terms in the series is only 3 like
  // (1^1000 + 2^1000 + 3^1000).  In this case, falling factorials in many
  // terms will be 0.  There is no point in computing the coefficients for
  // these terms.
  if (maxN < maxNumCoefficients) {
    maxNumCoefficients = maxN + 1;
  }

  for (k = 0; k < maxNumCoefficients; k++) {
    coeffs.push_back(0);
  }

  for (i = 0; i <= m; i++) {
    for (k = 0; k < maxNumCoefficients; k++) {
      if (i == k) {
        coeffs[k] = 1;
      } else if (i > 0 && k > 0) {
        // Calculate T(2*i, 2*k) = T(2*i-2,2*k-2) + k*k*T(2*i-2, 2*k)
        // Use the previous values in coeffs[k] for all the computation before
        // setting the new value.
        temp = coeffs[k]*k*k + T_2_2;
        T_2_2 = coeffs[k];
        coeffs[k] = temp;
      } else {
        coeffs[k] = 0;
      }
    }
    T_2_2 = coeffs[0];
  }
  return coeffs;
}

void CentralFactorialPowerSum :: printFallingFactorial(long start,
                                                long numTerms, ostream & out) {
  for(long i = 0; i < numTerms; i++) {
    if (i == 0) {
      out << "(n + " << start << ')';
    } else {
      long diff = start - i;
      if (diff == 0) {
        out << 'n';
      } else if (diff > 0) {
        out << "(n + " << diff << ')';
      } else {
        out << "(n - " << -diff << ')';
      }
    }
  }
}
