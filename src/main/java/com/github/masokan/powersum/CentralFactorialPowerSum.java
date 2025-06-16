package com.github.masokan.powersum;

import java.math.BigInteger;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class CentralFactorialPowerSum implements PowerSum {
  public CentralFactorialPowerSum() {
  }

  @Override
  public List<RationalNumber> getCoefficients(long power) {
    List<BigInteger> coeffs = getCoefficients(power, power);
    return PowerSum.toRationalNumbers(coeffs);
  }

  @Override
  public void printSumFormula(long power, PrintStream out) {
    if (power < 0) {
      return;
    }

    boolean evenPower = ((power & 1) == 0);

    if (power > 0) {
      List<BigInteger> coeffs = getCoefficients(power, power);
      // Coefficient 0 is always 0 - so we start the loop at index 1
      for (int k = 1; k < coeffs.size(); k++) {
        BigInteger coeff = coeffs.get(k);
        if (! BigInteger.ONE.equals(coeff)) {
          out.print(coeff);
        }
        if (evenPower) {
          out.print("(2n + 1)");
        }
        printFallingFactorial(k, 2*k, out);
        if (evenPower) {
          out.printf("/%d", 2*(2*k + 1));
        } else {
          out.printf("/%d", 2*k);
        }
        if (k != coeffs.size() - 1) {
          out.print(" + ");
        }
      }
    } else {
      // Special case for 0
      out.print("(n + 1)");
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

    if (power > 0) {
      before = PowerSum.getThreadCpuTime();
      List<BigInteger> coeffs = getCoefficients(power, n);
      after = PowerSum.getThreadCpuTime();
      stat[0] = after - before;
      before = PowerSum.getThreadCpuTime();
      boolean evenPower = ((power & 1) == 0);
      BigInteger fallingFactorial = BigInteger.ONE;

      // Coefficient 0 is always 0 - so we start the loop at index 1
      for (int k = 1; k < coeffs.size(); k++) {
        BigInteger coeff = coeffs.get(k);
        fallingFactorial = fallingFactorial.multiply(BigInteger.valueOf(n + k))
                           .multiply(BigInteger.valueOf(n - k + 1));
        if (evenPower) {
          sum = sum.add(coeff.multiply(fallingFactorial)
                             .multiply(BigInteger.valueOf(2*n + 1))
                             .divide(BigInteger.valueOf(2*(2*k + 1))));
        } else {
          sum = sum.add(coeff.multiply(fallingFactorial)
                             .divide(BigInteger.valueOf(2*k)));
        }
      }
    } else {
      stat[0] = 0;
      before = PowerSum.getThreadCpuTime();
      // Special case - not handled by the formula
      sum = BigInteger.valueOf(n + 1);
    }
    after = PowerSum.getThreadCpuTime();
    stat[1] = after - before;
    return sum;

  }

  /**
   * The Central Factorial Numbers of second kind are as defined below:
   * T(2m, 2k) = k*k*T(2m - 2, 2k) + T(2m - 2, 2k - 2)
   * T(2m, 2m) = 1
   * We need only even coefficients.
   */
  private List<BigInteger> getCoefficients(long power, long maxN) {
    ArrayList<BigInteger> coeffs = new ArrayList<>();
    if (power < 0) {
      return coeffs;
    }

    // Get half of power.  For odd power, round up half power
    long m = (power >> 1) + (power & 1);

    long maxNumCoefficients = m + 1; /* Array is 0 based, so add 1 */

    /* We do not compute coefficients more than what is needed.  Example: power
     * is 1000 and the number of terms in the series is only 3 like
     * (1^1000 + 2^1000 + 3^1000).  In this case, falling factorials in many
     * terms will be 0.  There is no point in computing the coefficients for
     * these terms.
     */
    if (maxN < maxNumCoefficients) {
      maxNumCoefficients = maxN + 1;
    }

    BigInteger temp;
    BigInteger T_2_2 = BigInteger.ZERO;
    long k;

    for (k = 0; k < maxNumCoefficients; k++) {
      coeffs.add(BigInteger.ZERO);
    }

    for (long i = 0; i <= m; i++) {
      for (k = 0; k < maxNumCoefficients; k++) {
        if (i == k) {
          coeffs.set((int)(k), BigInteger.ONE);
        } else if (i > 0 && k > 0) {
          temp = BigInteger.valueOf(k).multiply(BigInteger.valueOf(k))
                 .multiply(coeffs.get((int)(k)))
                 .add(T_2_2);
          T_2_2 = coeffs.get((int)(k));
          coeffs.set((int)(k), temp);
        } else {
          coeffs.set((int)(k), BigInteger.ZERO);
        }
      }
      T_2_2 = coeffs.get(0);
    }
    return coeffs;
  }

  private void printFallingFactorial(int start, int numTerms, PrintStream out) {
    for(int i = 0; i < numTerms; i++) {
      switch(i) {
        case 0:
          out.printf("(n + %d)", start);
          break;
        default:
          int diff = start - i;
          if (diff == 0) {
            out.print("n");
          } else if (diff > 0) {
            out.printf("(n + %d)", diff);
          } else {
            out.printf("(n - %d)", -diff);
          }
          break;
      }
    }
  }

}
