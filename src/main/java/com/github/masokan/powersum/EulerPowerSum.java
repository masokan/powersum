package com.github.masokan.powersum;

import java.math.BigInteger;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class EulerPowerSum implements PowerSum {
  public EulerPowerSum() {
  }

  @Override
  public void printSumFormula(long power, PrintStream out) {
    if (power < 0) {
      return;
    }

    List<BigInteger> coeffs = getCoefficients(power, power);
    BigInteger powerFactorial;

    // Compute (power + 1)!
    powerFactorial = BigInteger.ONE;
    for (long i = 2; i <= power + 1; i++) {
      powerFactorial = powerFactorial.multiply(BigInteger.valueOf(i));
    }

    out.print("{ ");
    for (int j = 0; j < coeffs.size(); j++) {
      if (coeffs.get(j).compareTo(BigInteger.ONE) == 1) {
        // Coefficient is greater than 1
        out.print(coeffs.get(j));
      }
      if (coeffs.get(j).compareTo(BigInteger.ZERO) == 1) {
        // Coefficient is greater than 0
        printTerm(j, power + 1, out);
      }
      if (j < (coeffs.size() - 2)) {
        out.print(" + ");
      }
    }
    out.print(" }");
    if (powerFactorial.compareTo(BigInteger.ONE) != 0) {
      out.print("/" + powerFactorial);
    }
    out.println();
  }

  @Override
  public BigInteger computeSumWithTimeStat(long power, long n, long[] stat) {
    long before;
    long after;
    BigInteger sum = BigInteger.ZERO;

    if (power < 0 || n < 0) {
      stat[0] = 0;
      stat[1] = 0;
      return sum;
    }

    before = PowerSum.getThreadCpuTime();
    List<BigInteger> coeffs = getCoefficients(power, n);
    after = PowerSum.getThreadCpuTime();
    stat[0] = after - before;
    before = PowerSum.getThreadCpuTime();
    BigInteger temp = BigInteger.ZERO;
    boolean firstTime = true;
    for (int k = 0; k < coeffs.size(); k++) {
      if (n + k >= power) {
        if (firstTime) {
          firstTime = false;
          temp = nCr(n + k + 1, power + 1);
        } else {
          temp = temp.multiply(BigInteger.valueOf(n + k + 1))
                     .divide(BigInteger.valueOf(n + k - power));
        }
        sum = sum.add(coeffs.get(k).multiply(temp));
      }
    }
    after = PowerSum.getThreadCpuTime();
    stat[1] = after - before;
    return sum;
  }

  @Override
  public List<RationalNumber> getCoefficients(long power) {
    List<BigInteger> coeffs = getCoefficients(power, power);
    return PowerSum.toRationalNumbers(coeffs);
  }

  /**
   * The Euler numbers of first kind are as defined below:
   * E(i, j) = 1 when j = 0
   * E(i, j) = (j + 1)E(n - 1, j) + (n - j)E(n - 1, j - 1)
   */
  private List<BigInteger> getCoefficients(long power,
                                           long maxNumCoefficients) {
    ArrayList<BigInteger> coeffs = new ArrayList<>();

    if (power < 0) {
      return coeffs;
    }

    BigInteger temp;
    BigInteger E_i_1_j_1 = BigInteger.ONE;
    long k;
    long j;

    coeffs.add(BigInteger.ONE);
    for (k = 1; k < power + 1; k++) {
      coeffs.add(BigInteger.ZERO);
    }

    for (long i = 1; i <= power; i++) {
      boolean oddPower = ((i & 1) == 1);
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
          E_i_1_j_1 = coeffs.get((int)(j));
          coeffs.set((int)(j), BigInteger.ONE);
        } else {
          E_i_1_j_1 = E_i_1_j_1.multiply(BigInteger.valueOf(i - j));
          temp = coeffs.get((int)(j)).multiply(BigInteger.valueOf(j + 1))
                                     .add(E_i_1_j_1);
          E_i_1_j_1 = coeffs.get((int)(j));
          coeffs.set((int)(j), temp);
        }
      }
      j = halfLimit;
      // Since coefficients are mirror reflection w.r.t. the central point,
      // we can reverse copy whatever we initialized so far
      if (oddPower) {
        for (k = 1; j + k < i; k++) {
          coeffs.set((int)(j + k), coeffs.get((int)(j - k)));
        }
      } else {
        for (k = 1; j + k < i; k++) {
          coeffs.set((int)(j + k), coeffs.get((int)(j - k + 1)));
        }
      }
      coeffs.set((int)i, BigInteger.ZERO);
    }
    return coeffs;
  }

  private void printTerm(long start, long numTerms, PrintStream out) {
    start++;
    for(long j = 0; j < numTerms; j++, start--) {
      if (start > 0) {
        out.print("(n + " + start + ")");
      } else if (start == 0) {
        out.print("n");
      } else {
        out.print("(n - " + (-start) + ")");
      }
    }
  }
}
