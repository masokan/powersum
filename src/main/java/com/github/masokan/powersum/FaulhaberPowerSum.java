package com.github.masokan.powersum;

import java.math.BigInteger;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class FaulhaberPowerSum implements PowerSum {
  public FaulhaberPowerSum() {
  }

  @Override
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
  public List<RationalNumber> getCoefficients(long power) {

    ArrayList<RationalNumber> coeffs = new ArrayList<>();

    if (power < 0) {
      return coeffs;
    }

    long nLimit;
    boolean oddPower = ((power & 1) == 1);

    if (oddPower) {
      nLimit = (power + 1)/2;
    } else {
      nLimit = power/2 + 1;
    }
    ArrayList<BigInteger> firstRow = new ArrayList<>();
    ArrayList<BigInteger> nextRow = new ArrayList<>();
    BigInteger scaleBy;
    BigInteger nextScaleBy;

    // Create the first row.  Each row is created with an augmented column entry
    // of 1
    if (oddPower) {
      scaleBy = createRowForOddPower(nLimit, nLimit, firstRow);
    } else {
      scaleBy = createRowForEvenPower(nLimit, nLimit, firstRow);
    }

    coeffs.add(new RationalNumber(firstRow.get(firstRow.size() - 1), scaleBy));
    long pivotIndex = 1;
    // Create subsequent rows and generate one coefficient at a time
    for (long i = nLimit - 1; i >= 1; i --) {
      // Reinitialize the augmented column
      firstRow.set(firstRow.size() - 1, BigInteger.ZERO);
      BigInteger pivot = firstRow.get((int)pivotIndex);
      if (pivot.equals(BigInteger.ZERO)) {
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
      scaleBy = scaleBy.multiply(nextScaleBy);

      // Get the value in the pivot column in the next row
      BigInteger pivotInNext = nextRow.get((int)pivotIndex);

      // Multiply rest of the columns that follow the pivot in the first row by
      // next pivot column and the next row by pivot and subtract from first
      for (long j = pivotIndex + 1; j < nLimit + 1; j++) {
        BigInteger temp = firstRow.get((int)j).multiply(pivotInNext);
        temp = temp.subtract(pivot.multiply(nextRow.get((int)j)));
        firstRow.set((int)j, temp);
      }
      coeffs.add(new RationalNumber(firstRow.get((int)nLimit), scaleBy)
                 .canonicalize());
      pivotIndex++;
    }

    return coeffs;
  }

  @Override
  public void printSumFormula(long power, PrintStream out) {
    if (power > 0) {
      List<RationalNumber> coeffs = getCoefficients(power);

      if ((power & 1) == 0) {
        out.print("(2n + 1)");
      }
      long exponent = (power + 1)/2;
      out.print("{");
      for (int i = 0; i < coeffs.size(); i++) {
        if (i != 0) {
          out.print(" + ");
        }
        out.print("(" + coeffs.get(i) + ")N");
        if (exponent != 1) {
          out.print("^" + exponent);
        }
        exponent--;
      }
      out.println("}/2");
      out.println("where N = n(n + 1)");
    } else {
      if (power == 0) {
        out.println("(n + 1)");
      }
    }
  }

  @Override
  public BigInteger computeSumWithTimeStat(long power, long n, long[] stat) {
    long before;
    long after;
    BigInteger result = BigInteger.ZERO;

    if (power < 0 || n < 0) {
      stat[0] = 0;
      stat[1] = 0;
      return result;
    }

    if (power > 0) {
      before = PowerSum.getThreadCpuTime();
      List<RationalNumber> coeffs = getCoefficients(power);
      after = PowerSum.getThreadCpuTime();
      stat[0] = after - before;
      before = PowerSum.getThreadCpuTime();
      if (n > 0) {
        RationalNumber sum = RationalNumber.ZERO;
        // Do the multiplication in BigInteger to avoid any overflow in long
        BigInteger N = BigInteger.valueOf(n);
        N = N.multiply(BigInteger.valueOf(n + 1));
        BigInteger NPow = N;

        // For odd powers, the first term starts with N^2.  There is no N term
        if ((power & 1) == 1) {
          NPow = NPow.multiply(N);
        }
        for (int i = coeffs.size() - 1; i >= 0 ; i--) {
          sum = sum.add(coeffs.get(i).multiply(NPow)).canonicalize();
          NPow = NPow.multiply(N);
        }
        if ((power & 1) == 0) {
          sum = sum.multiply(BigInteger.valueOf(2*n + 1)).canonicalize();
        }
        result = sum.getNumerator().divide(BigInteger.valueOf(2));
      }
    } else {
      stat[0] = 0;
      before = PowerSum.getThreadCpuTime();
      result = BigInteger.valueOf(n + 1);
    }
    after = PowerSum.getThreadCpuTime();
    stat[1] = after - before;
    return result;
  }

  /**
   * Reverse the row and return the first non-zero column value
   */
  private BigInteger postProcessRow(List<BigInteger> row) {
    Collections.reverse(row);
    row.add(BigInteger.ONE);
    BigInteger scaleBy = BigInteger.ONE;
    for (int j = 0; j < row.size(); j++) {
      if (! row.get(j).equals(BigInteger.ZERO)) {
        scaleBy = row.get(j);
        break;
      }
    }
    return scaleBy;
  }

  private BigInteger createRowForEvenPower(long nLimit, long rowNum,
                                           List<BigInteger> row) {
    long i;
    long numAdded = 0;

    BigInteger ncr1 = BigInteger.ZERO;
    BigInteger ncr2 = BigInteger.ZERO;

    for (i = 2*rowNum - 1; i >= 1; i -= 2) {
      if (i > rowNum) {
        row.add(BigInteger.ZERO);
      } else {
        // Calculate nCr only once and incrementally update it in the loop to
        // avoid redundant computation that will increase the time by an order
        // of magnitude
        if (ncr1.equals(BigInteger.ZERO)) {
          ncr1 = nCr(rowNum, i);
        } else {
          ncr1 = ncr1.multiply(BigInteger.valueOf((i + 2)*(i + 1)))
                 .divide(BigInteger.valueOf((rowNum - i - 1)*(rowNum - i)));
        }
        if (ncr2.equals(BigInteger.ZERO)) {
          ncr2 = nCr(rowNum - 1, i);
        } else {
          ncr2 = ncr2.multiply(BigInteger.valueOf((i + 2)*(i + 1)))
                 .divide(BigInteger.valueOf((rowNum - i - 2)*(rowNum - i - 1)));
        }
        row.add(ncr1.add(ncr2));
      }
      numAdded++;
    }
    while(numAdded < nLimit) {
      row.add(BigInteger.ZERO);
      numAdded++;
    }
    return postProcessRow(row);
  }

  private BigInteger createRowForOddPower(long nLimit, long rowNum,
                                          List<BigInteger> row) {
    long i;
    long numAdded = 0;

    BigInteger ncr = BigInteger.ZERO;

    for (i = 2*rowNum - 1; i >= 1; i -= 2) {
      if (i > rowNum) {
        row.add(BigInteger.ZERO);
      } else {
        // Calculate nCr only once and incrementally update it in the loop to
        // avoid redundant computation that will increase the time by an order
        // of magnitude
        if (ncr.equals(BigInteger.ZERO)) {
          ncr = nCr(rowNum, i);
        } else {
          ncr = ncr.multiply(BigInteger.valueOf((i + 2)*(i + 1)))
                .divide(BigInteger.valueOf((rowNum - i - 1)*(rowNum - i)));
        }
        row.add(ncr);
      }
      numAdded++;
    }
    while(numAdded < nLimit) {
      row.add(BigInteger.ZERO);
      numAdded++;
    }
    return postProcessRow(row);
  }


}
