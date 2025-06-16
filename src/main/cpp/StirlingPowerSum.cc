#include <iostream>

using std::endl;

#include "StirlingPowerSum.h"

StirlingPowerSum :: StirlingPowerSum()
                  : PowerSum() {
}

vector<mpq_class> StirlingPowerSum :: getCoefficients(long power) {
  vector<mpq_class> outCoeffs;
  vector<mpz_class> coeffs = getCoefficients(power, power + 1);
  for (size_t i = 0; i < coeffs.size(); i++) {
    outCoeffs.push_back(mpq_class(coeffs[i]));
  }
  return outCoeffs;
}

void StirlingPowerSum :: printSumFormula(long power, ostream &out) {
  if (power < 0) {
    return;
  }

  vector<mpz_class> coeffs = getCoefficients(power, power + 1);
  bool firstTime = true;
  for (size_t t = 0; t < coeffs.size(); t++) {
    mpz_class coeff = coeffs[t];
    if (coeff != 0) {
      if (firstTime) {
        out << "   ";
        firstTime = false;
      } else {
        out << " + ";
      }
      if (coeff > 1) {
        // Coefficient is greater than 1
        out << coeff;
      }
      printFactors(t, out);
      if (t > 0) {
        out << '/' << (t + 1);
      }
    }
  }
  out << endl;
}

mpz_class StirlingPowerSum :: computeSumWithTimeStat(long power, long n, 
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
  vector<mpz_class> coeffs = getCoefficients(power, n + 1);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  mpz_class fallingFactorial = n + 1;
  mpz_class factor = 0;
  long numTermsToCompute = coeffs.size();
  if (numTermsToCompute > (n + 1)) {
    numTermsToCompute = n + 1;
  }
  for (long t = 0; t < numTermsToCompute; t++) {
    // We know that the following division is exact and we do this first
    // before the multiplication to avoid generating a large intermediate
    // value
    factor = fallingFactorial/(t + 1);
    sum += coeffs[t]*factor;;
    fallingFactorial *= (n - t);
  }
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));
  return sum;
}

StirlingPowerSum :: ~StirlingPowerSum() {
}

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
vector<mpz_class> StirlingPowerSum :: getCoefficients(long power,
                                              long maxNumCoefficients) {
  vector<mpz_class> coeffs;

  if (power < 0) {
    return coeffs;
  }

  if (power > 0) {
    maxNumCoefficients--;
    // Initialize all coefficients to 0
    for (long i = 0; i <= power; i++) {
      coeffs.push_back(0);
    }

    // Anchors to hold values from previous iteration
    mpz_class Smj = 0;
    mpz_class Sm_1j_1 = 1; // S(0, 0)
  
    for (long currentPower = 1; currentPower <= power; currentPower++) {
      long termMax
       = (currentPower > maxNumCoefficients)? maxNumCoefficients : currentPower;
      // Initialize other coefficients for the current power
      for(long term = 1; term <= termMax; term++) {
        Smj = term*coeffs[term] + Sm_1j_1;
        Sm_1j_1 = coeffs[term];
        coeffs[term] = Smj;
      }
      Sm_1j_1 = coeffs[0];
    }
  } else {
    coeffs.push_back(1);
  }
  return coeffs;
}

void StirlingPowerSum :: printFactors(long term, ostream & out) {
  for(long i = 0; i <= term; i++) {
    switch(i) {
      case 0:
        out << "(n + 1)";
        break;
      case 1:
        out << 'n';
        break;
      default:
        out << "(n - "  << (i - 1) << ')';
        break;
    }
  }
}
