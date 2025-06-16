#include <iostream>

using std::endl;

#include "FaulhaberPowerSum.h"

FaulhaberPowerSum :: FaulhaberPowerSum()
                  : PowerSum() {
}

/**
 * There is no simple recurrence relation to generate the coefficients.
 * According to A. W. F. Edwards
 * (http://www.pietrocola.eu/Fontecchio2019/A%20quick%20route%20to%20sums%20of%20powers%20by%20A.W.F.Edwards%20(1).pdf),
 * they can be obtained by matrix inversion where matrix rows are initialized
 * using the methods outlined in the paper.  In general, a matrix requires
 * O(m^2) storage and the inversion requires O(m^p) time where p > 2.
 * Instead of matrix inversion, we consider the problem as solving a system of
 * linear equtaions.  So, we contruct an augmented matrix and transform it to
 * reduced row echelon (rre) form.  In our case, we notice that if we properly
 * construct the matrix, it can be a row echelon matrix.
 * All we need to do is transform it to a rre matrix and solve for the
 * coefficients.  The transformation can be achieved sequentially one row at a
 * time.  There is no need to hold the entire matrix in memory thus reducing
 * the storage requirement to O(m).  The time complexity is also brought down
 * to O(m^2).
 */
vector<mpq_class> FaulhaberPowerSum :: getCoefficients(long power) {
  vector<mpq_class> coefficients;

  if (power < 0) {
    return coefficients;
  }

  vector<mpz_class> firstRow;
  vector<mpz_class> nextRow;
  mpz_class scaleBy;
  mpz_class nextScaleBy;

  long nLimit;
  bool oddPower = false;
  if (power & 1) {
    oddPower = true;
    nLimit = (power + 1)/2;
  } else {
    nLimit = power/2 + 1;
  }
  // Create the first row.  Each row is created with an augmented column entry
  // of 1
  if (oddPower) {
    scaleBy = createRowForOddPower(nLimit, nLimit, firstRow);
  } else {
    scaleBy = createRowForEvenPower(nLimit, nLimit, firstRow);
  }

  coefficients.push_back(mpq_class(firstRow.back(), scaleBy));
  long pivotIndex = 1;
  // Create subsequent rows and generate one coefficient at a time
  for (long i = nLimit - 1; i >= 1; i --) {
    // Reinitialize the augmented column
    firstRow.back() = 0;
    mpz_class pivotInFirst = firstRow[pivotIndex];
    if (pivotInFirst == 0) {
      // Once we reach a 0, all other columns will be zero.  This means all
      // other coefficients will be zero
      break;
    }
    nextRow.clear();
    if (oddPower) {
      nextScaleBy = createRowForOddPower(nLimit, i, nextRow);
    } else {
      nextScaleBy = createRowForEvenPower(nLimit, i, nextRow);
    }
    scaleBy *= nextScaleBy;

    // Get the value in the pivot column in the next row
    mpz_class pivotInNext = nextRow[pivotIndex];
 
    // Multiply rest of the columns that follow the pivot in the first row by
    // next pivot column and the next row by pivot and subtract from first
    for (size_t j = pivotIndex + 1; j < nextRow.size(); j++) {
      firstRow[j] *= pivotInNext;
      firstRow[j] -= pivotInFirst*nextRow[j];
    }
    mpq_class coeff = mpq_class(firstRow.back(), scaleBy);
    coeff.canonicalize();
    coefficients.push_back(coeff);
    pivotIndex++;
  }

  return coefficients;
}

void FaulhaberPowerSum :: printSumFormula(long power, ostream &out) {
  if (power > 0) {
    vector<mpq_class> coeffs = getCoefficients(power);
    if ((power & 1) == 0) {
      out << "(2n + 1)";
    }
    long exponent = (power + 1)/2;
    out << '{';
    for (size_t i = 0; i < coeffs.size(); i++) {
      if (i != 0) {
        out << " + ";
      }
      out << '(' << coeffs[i] << ")N";
      if (exponent != 1) {
        out << '^' << exponent;
      }
      exponent--;
    }
    out << "}/2" << endl;
    out << "where N = n(n + 1)" << endl;
  } else {
    if (power == 0) {
      out << "(n + 1)" << endl;
    }
  }
}

mpz_class FaulhaberPowerSum :: computeSumWithTimeStat(long power, long n, 
                                                      vector<long> & stat) {
  struct timespec before;
  struct timespec after;

  stat.clear();
  mpq_class sum = 0;

  if (power < 0 || n < 0) {
    stat.push_back(0);
    stat.push_back(0);
    return sum.get_num();
  }

  if (power > 0) {
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    vector<mpq_class> coeffs = getCoefficients(power);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
    stat.push_back(computeCpuTime(before, after));
    if (n > 0) {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
      mpz_class N = n;
      N *= (n + 1); // To avoid overflow in long, do the multiplication in mpz
      mpz_class NPow = N;
      // For odd powers, the first term starts with N^2.  There is no N term
      if ((power & 1) == 1) {
        NPow *= N;
      }
      long coeffSize = (long)coeffs.size();
      for (long i = coeffSize - 1; i >= 0; i--) {
        sum += coeffs[i]*NPow;
        NPow *= N;
      }
      sum /= 2;
      if ((power & 1) == 0) {
        sum *= (2*n + 1);
      }

    }
      
  } else {
    stat.push_back(0);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    sum = n + 1;
  }
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  stat.push_back(computeCpuTime(before, after));
  return sum.get_num();
}

FaulhaberPowerSum :: ~FaulhaberPowerSum() {
}

/**
 * Reverse the row and return the first non-zero column value
 */
mpz_class FaulhaberPowerSum :: postProcessRow(vector<mpz_class> & row) {
  std::reverse(row.begin(), row.end());
  // Augment the row with 1 which corresponds to the coefficient of the
  // power of N (Note: N = n(n+1))
  row.push_back(1);
  mpz_class scaleBy = 1;
  for (size_t j = 0; j < row.size(); j++) {
    if (row[j] != 0) {
      scaleBy = row[j];
      break;
    }
  }
  return scaleBy;
}

mpz_class FaulhaberPowerSum :: createRowForEvenPower(long nLimit, long rowNum,
                                                     vector<mpz_class> & row) {
  long numAdded = 0;
  mpz_class ncr1 = 0;
  mpz_class ncr2 = 0;
  for (long i = 2*rowNum - 1; i >= 1; i -= 2) {
    if (i > rowNum) {
      row.push_back(0);
    } else {
      // Calculate nCr only once and incrementally update it in the loop to
      // avoid redundant computation that will increase the time by an order
      // of magnitude
      if (ncr1 == 0) {
        ncr1 = nCr(rowNum, i);
      } else {
        ncr1 *= (i + 2)*(i + 1);
        ncr1 /= (rowNum - i - 1)*(rowNum - i);
      }
      if (ncr2 == 0) {
        ncr2 = nCr(rowNum - 1, i);
      } else {
        ncr2 *= (i + 2)*(i + 1);
        ncr2 /= (rowNum - i - 2)*(rowNum - i - 1);
      }
      row.push_back(ncr1 + ncr2);
    }
    numAdded++;
  }
  while(numAdded < nLimit) {
    row.push_back(0);
    numAdded++;
  }
  return postProcessRow(row);
}

mpz_class FaulhaberPowerSum :: createRowForOddPower(long nLimit, long rowNum,
                                                    vector<mpz_class> & row) {

  long numAdded = 0;
  mpz_class ncr = 0;
  for (long i = 2*rowNum - 1; i >= 1; i -= 2) {
    if (i > rowNum) {
      row.push_back(0);
    } else {
      // Calculate nCr only once and incrementally update it in the loop to
      // avoid redundant computation that will increase the time by an order
      // of magnitude
      if (ncr == 0) {
        ncr = nCr(rowNum, i);
      } else {
        ncr *= (i + 2)*(i + 1);
        ncr /= (rowNum - i - 1)*(rowNum - i);
      }
      row.push_back(ncr);
    }
    numAdded++;
  }
  while(numAdded < nLimit) {
    row.push_back(0);
    numAdded++;
  }
  return postProcessRow(row);
}
