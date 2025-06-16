package com.github.masokan.powersum;

import java.math.BigInteger;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class BernoulliPowerSum implements PowerSum {
  public BernoulliPowerSum() {
  }

  @Override
  /* Bernoulli numbers are defined by the following recurrence relation:
   * B(0) = 1
   * B(m) = -(Binom((m + 1), 0)B(0) + Binom((m + 1), 1)B(1)
   *          + ... + Binom((m + 1), (m - 1))B(m - 1))
   * where Binom(i, j) is the binomial coefficient which evaluates to:
   * i!/{(i - j)!j!}
   */
  public List<RationalNumber> getCoefficients(long power) {
    ArrayList<RationalNumber> coeffs = new ArrayList<>();

    if (power < 0) {
      return coeffs;
    }

    long i;
    // Initialize B(0)
    coeffs.add(RationalNumber.ONE);
    if (power > 0) {
      // Initialize B(1)
      coeffs.add(new RationalNumber(-1, 2));

      for (i = 2; i <= power; i++) {
        if ((i & 1) == 1) {
          // Odd coefficients above 1 are 0s
          coeffs.add(RationalNumber.ZERO);
        } else {
          coeffs.add(computeNextCoefficient(coeffs, i));
        }
      }
    }
    return coeffs;
  }

  @Override
  public void printSumFormula(long power, PrintStream out) {
    if (power < 0) {
      return;
    }

    List<RationalNumber> coeffs = getCoefficients(power);
    BigInteger binom = BigInteger.ONE;
    long binomN = power + 1;
    long binomR = 1;
    long p = power + 1;

    out.print("{ ");
    for (int i = 0; i < coeffs.size(); i++) {
      if ((i & 1) == 0 || i == 1) {
        if (i != 0) {
          out.print(" + ");
        }
        if (! RationalNumber.ONE.equals(coeffs.get(i))) {
          out.print("(" + coeffs.get(i) + ")");
        }
        if (! BigInteger.ONE.equals(binom)) {
          out.print("" + binom);
        }
        out.print("(n + 1)");
        if (p != 1) {
          out.print("^" + p);
        }
      }
      p--;
      binom = binom.multiply(BigInteger.valueOf(binomN))
                   .divide(BigInteger.valueOf(binomR));
      binomN--;
      binomR++;
    }
    
    out.print(" }");
    if (power > 0) {
      out.print("/" + (power + 1));
    }
    out.println();
  }

  @Override
  public BigInteger computeSumWithTimeStat(long power, long n, long[] stat) {
    long before;
    long after;
    RationalNumber sum = RationalNumber.ZERO;

    if (power < 0 || n < 0) {
      stat[0] = 0;
      stat[1] = 0;
      return sum.getNumerator();
    }

    before = PowerSum.getThreadCpuTime();
    List<RationalNumber> coeffs = getCoefficients(power);
    after = PowerSum.getThreadCpuTime();
    stat[0] = after - before;
    before = PowerSum.getThreadCpuTime();
    BigInteger bign = BigInteger.valueOf(n + 1);
    BigInteger pow = bign;
    long binomN = power + 1;
    long binomR = 1;
    BigInteger binom = BigInteger.valueOf(binomN);

    // We compute the terms in reverse order to avoid an additional division
    // operation in the loop
    for (int i = coeffs.size() - 1; i >= 0; i--) {
      if ((i & 1) == 0 || i == 1) {
        sum = sum.add(coeffs.get(i).multiply(binom).multiply(pow)
                            .canonicalize());
      }
      pow = pow.multiply(bign);
      binomN--;
      binomR++;
      binom = binom.multiply(BigInteger.valueOf(binomN))
                   .divide(BigInteger.valueOf(binomR));
    }
    sum = sum.divide(BigInteger.valueOf(power + 1)).canonicalize();
    after = PowerSum.getThreadCpuTime();
    stat[1] = after - before;
    return sum.getNumerator();
  }

  private RationalNumber computeNextCoefficient(List<RationalNumber> current,
                                                long m) {
    BigInteger binom = BigInteger.ONE;
    RationalNumber coeff = RationalNumber.ZERO;

    for (long k = 0; k < m; k++) {
      if (((k & 1) == 0) || k == 1) {
        // Even coefficient or the first odd one
        coeff = coeff.add(current.get((int)(k)).multiply(binom));
      }
      binom = binom.multiply(BigInteger.valueOf((m + 1 - k)));
      binom = binom.divide(BigInteger.valueOf((k + 1)));
    }
    coeff = coeff.negate().divide(binom).canonicalize();
    return coeff;

  }

}
