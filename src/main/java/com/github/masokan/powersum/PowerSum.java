package com.github.masokan.powersum;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.math.BigInteger;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public interface PowerSum {
  ThreadMXBean threadMXBean = ManagementFactory.getThreadMXBean();

  /**
   * To get the coefficients in the sum formula.
   * @param power - desired power
   * @return a list of coefficients
   */
  public List<RationalNumber> getCoefficients(long power);

  /**
   * To print the sum formula for a given power
   * @param power - desired power
   * @param out - output stream where the formula should be printed
   */
  public void printSumFormula(long power, PrintStream out);

  /**
   * To compute the sum for a specific power and the number of terms and obtain
   * CPU times
   * @param power - desired power
   * @param numTerms - number of terms
   * @param stat - array in which  the coefficient initialization time and
   *               summation time (in nanosec.) are returned
   * @return sum computed
   */
  public BigInteger computeSumWithTimeStat(long power, long n, long[] stat);

  /**
   * To compute the sum for a specific power and the number of terms
   * @param power - desired power
   * @param numTerms - number of terms
   * @return sum computed
   */
  public default BigInteger computeSum(long power, long n) {
    long[] stat = new long[2];
    return computeSumWithTimeStat(power, n, stat);
  }

  /**
   * To compute the sum for a specific power and the number of terms using the
   * simple implementation (series summation)
   * @param power - desired power
   * @param numTerms - number of terms
   * @return sum computed
   */
  public default BigInteger computeSumUsingSeries(long power, long n) {
    BigInteger sum = BigInteger.ZERO;
    if (power < 0 || n < 0) {
      return sum;
    }
    for (long i = 0; i <= n; i++) {
      BigInteger base = BigInteger.valueOf(i);
      sum = sum.add(base.pow((int)power));
    }
    return sum;
  }

  public default BigInteger nCr(long n, long r) {
    long num = n;
    long i;

    if (r > n) {
      return BigInteger.ZERO;
    }
    if (r > n/2) {
      r = n - r;
    }
    if (r == 0) {
      return BigInteger.ONE;
    }
    BigInteger result = BigInteger.valueOf(num);

    for (i = 2; i <= r; i++) {
      num--;
      result = result.multiply(BigInteger.valueOf(num))
              .divide(BigInteger.valueOf(i));
    }
    return result;
  }

  static public List<RationalNumber>
      toRationalNumbers(List<BigInteger> coeffs) {
    List<RationalNumber> outCoeffs = new ArrayList<>();
    for (int i = 0; i < coeffs.size(); i++) {
      outCoeffs.add(new RationalNumber(coeffs.get(i)));
    }
    return outCoeffs;
  }

  static public long getThreadCpuTime() {
    return threadMXBean.getThreadCpuTime(Thread.currentThread().getId());
  }
}
