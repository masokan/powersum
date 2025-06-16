#include <stdlib.h>

#include "faulhaber_power_sum.h"

/**
 * Reverse the row and scale it by the first non-zero column so that the
 * leading column is 1.  This means the matrix will be in row echelon format
 */
static void postProcessRow(mpz_t * row, long rowSize, mpz_t scaleBy) {
  /* Reverse the row up to rowSize - 1 elements.  The last element will be
     initialized to 1 after reversing */
  long left = 0;
  long right = rowSize - 2;
  while (left < right) {
    /* Exchange the elements */
    mpz_swap(row[left], row[right]);
    /* Adjust left and right */
    left++;
    right--;
  }

  /* Augment the row with 1 which corresponds to the coefficient of the
     power of N (Note: N = n(n+1)) */
  mpz_set_ui(row[rowSize - 1], 1UL);
  mpz_set_ui(scaleBy, 1UL);
  for (long j = 0; j < rowSize; j++) {
    if (mpz_cmp_ui(row[j], 0) != 0) {
      mpz_set(scaleBy, row[j]);
      break;
    }
  }
}

static void createRowForEvenPower(long nLimit, long rowNum, mpz_t * row,
                                  mpz_t scaleBy) {
  long i;
  long numAdded = 0;

  mpz_t ncr1;
  mpz_init_set_ui(ncr1, 0UL);
  mpz_t ncr2;
  mpz_init_set_ui(ncr2, 0UL);

  mpz_t tempz;
  mpz_init(tempz);

  for (i = 2*rowNum - 1; i >= 1; i -= 2) {
    if (i > rowNum) {
      mpz_set_ui(row[numAdded], 0UL);
    } else {
      /* Calculate nCr only once and incrementally update it in the loop to
         avoid redundant computation that will increase the time by an order
         of magnitude */
      if (mpz_cmp_ui(ncr1,  0UL) == 0) {
        nCr(rowNum, i, ncr1);
      } else {
        mpz_mul_si(ncr1, ncr1, (i + 2)*(i + 1));
        mpz_divexact_ui(ncr1, ncr1, (rowNum - i - 1)*(rowNum - i));
      }
      if (mpz_cmp_ui(ncr2, 0UL) == 0) {
        nCr(rowNum - 1, i, ncr2);
      } else {
        mpz_mul_si(ncr2, ncr2, (i + 2)*(i + 1));
        mpz_divexact_ui(ncr2, ncr2, (rowNum - i - 2)*(rowNum - i - 1));
      }
      mpz_add(tempz, ncr1, ncr2);
      mpz_set(row[numAdded], tempz);
    }
    numAdded++;
  }
  while(numAdded < nLimit) {
    mpz_set_ui(row[numAdded], 0UL);
    numAdded++;
  }
  postProcessRow(row, nLimit + 1, scaleBy);
  mpz_clear(tempz);
  mpz_clear(ncr1);
  mpz_clear(ncr2);
}

static void createRowForOddPower(long nLimit, long rowNum, mpz_t * row,
                                 mpz_t scaleBy) {
  long numAdded = 0;

  mpz_t ncr;
  mpz_init_set_ui(ncr, 0UL);

  for (long i = 2*rowNum - 1; i >= 1; i -= 2) {
    if (i > rowNum) {
      mpz_set_ui(row[numAdded], 0UL);
    } else {
      /* Calculate nCr only once and incrementally update it in the loop to
         avoid redundant computation that will increase the time by an order
         of magnitude */
      if (mpz_cmp_ui(ncr,  0UL) == 0) {
        nCr(rowNum, i, ncr);
      } else {
        mpz_mul_si(ncr, ncr, (i + 2)*(i + 1));
        mpz_divexact_ui(ncr, ncr, (rowNum - i - 1)*(rowNum - i));
      }
      mpz_set(row[numAdded], ncr);
    }
    numAdded++;
  }
  while(numAdded < nLimit) {
    mpz_set_ui(row[numAdded], 0UL);
    numAdded++;
  }
  postProcessRow(row, nLimit + 1, scaleBy);
  mpz_clear(ncr);
}

static void freeRow(mpz_t * row, long size) {
  for (long i = 0; i < size; i++) {
    mpz_clear(row[i]);
  }
  free(row);
}

void freeFaulhaberCoeffs(mpq_t * coeffs, long size) {
  for (long i = 0; i < size; i++) {
    mpq_clear(coeffs[i]);
  }
  free(coeffs);
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
void getFaulhaberCoefficients(long power, mpq_t ** coeffList,
                              long * numCoeffs) {
  long nLimit;
  unsigned long oddPower = 0;
  mpz_t scaleBy;
  mpz_t nextScaleBy;
  mpz_t pivot;
  mpz_t tempz;
  mpq_t tempq;
  mpq_t scaleq;
  long j;

  if (power < 0) {
    *coeffList = NULL;
    *numCoeffs = 0;
    return;
  }

  mpz_init(scaleBy);
  mpz_init(nextScaleBy);
  mpz_init(pivot);
  mpz_init(tempz);
  mpq_init(tempq);
  mpq_init(scaleq);

  if (power & 1) {
    oddPower = 1;
    nLimit = (power + 1)/2;
  } else {
    nLimit = power/2 + 1;
  }
  mpq_t * coeffs = (mpq_t *)malloc(nLimit*sizeof(mpq_t));
  mpz_t * firstRow = (mpz_t *)malloc((nLimit + 1)*sizeof(mpz_t));
  mpz_t * nextRow = (mpz_t *)malloc((nLimit + 1)*sizeof(mpz_t));

  for (j = 0; j < nLimit; j++) {
    mpq_init(coeffs[j]);
    mpz_init(firstRow[j]);
    mpz_init(nextRow[j]);
  }

  /* Initialize augmented column */
  mpz_init(firstRow[nLimit]);
  mpz_init(nextRow[nLimit]);

  /* Create the first row.  Each row is created with an augmented column entry
   * of 1 */

  if (oddPower) {
    createRowForOddPower(nLimit, nLimit, firstRow, scaleBy);
  } else {
    createRowForEvenPower(nLimit, nLimit, firstRow, scaleBy);
  }

  mpq_set_z(scaleq, scaleBy);
  mpq_set_z(coeffs[0], firstRow[nLimit]);
  mpq_div(coeffs[0], coeffs[0], scaleq);

  long pivotIndex = 1;
  /* Create subsequent rows and generate one coefficient at a time */
  for (long i = nLimit - 1; i >= 1; i --) {
    /* Reinitialize the augmented column */
    mpz_set_si(firstRow[nLimit], 0);
    mpz_set(pivot, firstRow[pivotIndex]);
    if (mpz_cmp_ui(pivot, 0) == 0) {
      /* Once we reach a 0, all other columns will be zero.  This means all
       * other coefficients will be zero */
      break;
    }
    if (oddPower) {
      createRowForOddPower(nLimit, i, nextRow, nextScaleBy);
    } else {
      createRowForEvenPower(nLimit, i, nextRow, nextScaleBy);
    }
    mpz_mul(scaleBy, scaleBy, nextScaleBy);
    mpq_set_z(scaleq, scaleBy);

    /* Multiply rest of the columns that follow the pivot in the first row by
       next pivot column and the next row by pivot and subtract from first */
    for (j = pivotIndex + 1; j < nLimit + 1; j++) {
      mpz_mul(firstRow[j], firstRow[j], nextRow[pivotIndex]);
      mpz_mul(tempz, pivot, nextRow[j]);
      mpz_sub(firstRow[j], firstRow[j], tempz);
    }
    mpq_set_z(tempq, firstRow[nLimit]);
    mpq_div(tempq, tempq, scaleq);
    mpq_set(coeffs[pivotIndex], tempq);
    pivotIndex++;
  }

  mpq_t * finalCoeffs = (mpq_t *)malloc(pivotIndex*sizeof(mpq_t));
  for (j = 0; j < pivotIndex; j++) {
    mpq_init(finalCoeffs[j]);
    mpq_set(finalCoeffs[j], coeffs[j]);
  }

  /* Release storages */
  freeRow(firstRow, nLimit + 1);
  freeRow(nextRow, nLimit + 1);
  freeFaulhaberCoeffs(coeffs, nLimit);

  mpq_clear(scaleq);
  mpq_clear(tempq);
  mpz_clear(tempz);
  mpz_clear(pivot);
  mpz_clear(nextScaleBy);
  mpz_clear(scaleBy);

  *coeffList = finalCoeffs;
  *numCoeffs = pivotIndex;
}

void printFaulhaberSumFormula(long power, FILE * out) {
  if (power > 0) {
    mpq_t * coeffs;
    long numCoeffs;
    mpz_t zero;
    mpz_init_set_ui(zero, 0UL);
    getFaulhaberCoefficients(power, &coeffs, &numCoeffs);

    if ((power & 1) == 0) {
      fprintf(out, "(2n + 1)");
    }
    long exponent = (power + 1)/2;
    fprintf(out, "{");
    for (long i = 0; i < numCoeffs && mpq_cmp_z(coeffs[i], zero) != 0; i++) {
      if (i != 0) {
        fprintf(out, " + ");
      }
      gmp_fprintf(out, "(%Qd)N", coeffs[i]);
      if (exponent != 1) {
        fprintf(out, "^%ld", exponent);
      }
      exponent--;
    }
    fprintf(out, "}/2\n");
    fprintf(out, "where N = n(n + 1)\n");
    mpz_clear(zero);
    freeFaulhaberCoeffs(coeffs, numCoeffs);
  } else {
    if (power == 0) {
      fprintf(out, "(n + 1)\n");
    }
  }
}

void computeFaulhaberSumWithTimeStat(long power, long n, mpz_t result,
                                     long * coeffInitTime, long * sumTime) {
  struct timespec before;
  struct timespec after;
  mpz_set_ui(result, 0UL);

  if (power < 0 || n < 0) {
    *coeffInitTime = 0;
    *sumTime = 0;
    return;
  }

  if (power > 0) {
    mpq_t * coeffs;
    long numCoeffs;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    getFaulhaberCoefficients(power, &coeffs, &numCoeffs);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
    *coeffInitTime = computeCpuTime(&before, &after);
    if (n > 0) {
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
      mpq_t sum;
      mpq_init(sum);
      mpq_set_ui(sum, 0UL, 1UL);

      mpz_t N;
      mpz_init_set_ui(N, n);
      mpz_mul_ui(N, N, n + 1); /* Do this multiplication in mpz to avoid
                                  overflow when we do n*(n+1) in long */
      mpz_t NPow;
      mpz_init(NPow);

      mpz_set(NPow, N);
      mpq_t tempq;
      mpq_init(tempq);

      /* For odd powers, the first term starts with N^2.  There is no N term */
      if ((power & 1) == 1) {
        mpz_mul(NPow, NPow, N);
      }
      for (long i = numCoeffs - 1; i >= 0; i--) {
        mpq_set_z(tempq, NPow);
        mpq_mul(tempq, tempq, coeffs[i]);
        mpq_add(sum, sum, tempq);
        mpz_mul(NPow, NPow, N);
      }
      if ((power & 1) == 0) {
        mpq_set_ui(tempq, (2*n + 1), 1);
        mpq_mul(sum, sum, tempq);
      }
      mpq_get_num(result, sum);
      mpz_divexact_ui(result, result, 2UL);
      freeFaulhaberCoeffs(coeffs, numCoeffs);
      mpq_clear(tempq);
      mpz_clear(NPow);
      mpz_clear(N);
      mpq_clear(sum);
    }
  } else {
    *coeffInitTime = 0;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
    mpz_set_ui(result, n + 1);
  }
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  *sumTime = computeCpuTime(&before, &after);
}

void computeFaulhaberSum(long power, long n, mpz_t result) {
  long cTime;
  long sTime;

  computeFaulhaberSumWithTimeStat(power, n, result, &cTime, &sTime);
}

void initFaulhaberPowerSum(PowerSum * psPtr) {
  initPowerSum(psPtr);
  psPtr->getCoefficients = getFaulhaberCoefficients;
  psPtr->printSumFormula = printFaulhaberSumFormula;
  psPtr->computeSumWithTimeStat = computeFaulhaberSumWithTimeStat;
  psPtr->computeSum = computeFaulhaberSum;
  psPtr->freeCoefficients = freeFaulhaberCoeffs;
}
