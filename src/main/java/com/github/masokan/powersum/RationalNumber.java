package com.github.masokan.powersum;

import java.math.BigInteger;

/**
 * A lean implementation of RationalNumber with enough methods for use within
 * power sum computation.
 */
public class RationalNumber implements Comparable<RationalNumber> {
  private BigInteger numerator;
  private BigInteger denominator;

  static final public RationalNumber ZERO = new RationalNumber(0);
  static final public RationalNumber ONE = new RationalNumber(1);

  public RationalNumber(BigInteger num, BigInteger denom) {
    if (denom.signum() == 0) {
      throw new IllegalArgumentException("Denominator cannot be zero");
    }
    numerator = num;
    denominator = denom;
    adjustSigns();
  }

  public RationalNumber(BigInteger num) {
    numerator = num;
    denominator = BigInteger.ONE;
  }

  public RationalNumber(int num) {
    this(BigInteger.valueOf(num));
  }

  public RationalNumber(int num, int denom) {
    this(BigInteger.valueOf(num), BigInteger.valueOf(denom));
  }

  public BigInteger getNumerator() {
    return numerator;
  }

  public BigInteger getDenominator() {
    return denominator;
  }

  public RationalNumber negate() {
    return new RationalNumber(numerator.negate(), denominator);
  }

  public RationalNumber add(BigInteger a) {
    BigInteger num = numerator.add(denominator.multiply(a));
    return new RationalNumber(num, denominator);
  }

  public RationalNumber add(RationalNumber a) {
    BigInteger num = numerator.multiply(a.denominator)
                              .add(denominator.multiply(a.numerator));
    BigInteger denom = denominator.multiply(a.denominator);
    return new RationalNumber(num, denom);
  }

  public RationalNumber subtract(BigInteger s) {
    BigInteger num = numerator.subtract(denominator.multiply(s));
    return new RationalNumber(num, denominator);
  }

  public RationalNumber subtract(RationalNumber s) {
    BigInteger num = numerator.multiply(s.denominator)
                              .subtract(denominator.multiply(s.numerator));
    BigInteger denom = denominator.multiply(s.denominator);
    return new RationalNumber(num, denom);
  }

  public RationalNumber multiply(BigInteger m) {
    BigInteger num = numerator.multiply(m);
    return new RationalNumber(num, denominator);
  }

  public RationalNumber multiply(RationalNumber m) {
    BigInteger num = numerator.multiply(m.numerator);
    BigInteger denom = denominator.multiply(m.denominator);
    return new RationalNumber(num, denom);
  }

  public RationalNumber divide(BigInteger m) {
    if (m.signum() == 0) {
      throw new IllegalArgumentException("Divide by zero");
    }
    BigInteger denom = denominator.multiply(m);
    return new RationalNumber(numerator, denom);
  }

  public RationalNumber divide(RationalNumber m) {
    if (m.numerator.signum() == 0) {
      throw new IllegalArgumentException("Divide by zero");
    }
    BigInteger num = numerator.multiply(m.denominator);
    BigInteger denom = denominator.multiply(m.numerator);
    return new RationalNumber(num, denom);
  }

  @Override
  public boolean equals(Object other) {
    if (this == other) {
      return true;
    }
    boolean result = false;
    if (other != null && other instanceof RationalNumber) {
      RationalNumber o = (RationalNumber)other;
      int thisSign = getSign();
      int otherSign = o.getSign();
      if (thisSign == otherSign) {
        BigInteger m1 = numerator.multiply(o.denominator);
        BigInteger m2 = denominator.multiply(o.numerator);
        result = m1.equals(m2);
      }
    }
    return result;
  }

  public int getSign() {
    return numerator.signum();
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append(numerator.toString());
    if (denominator.compareTo(BigInteger.ONE) != 0) {
      sb.append("/").append(denominator.toString());
    }
    return sb.toString();
  }

  @Override
  public int compareTo(RationalNumber other) {
    int thisSign = getSign();
    int otherSign = other.getSign();
    if (thisSign != otherSign) {
      if (thisSign == -1) {
        return -1;
      } else {
        return 1;
      }
    } else if (thisSign == 0) {
      return 0;
    } else {
      return this.subtract(other).getSign();
    }
  }

  /** Reduce this rational number by dividing the numerator and denominator
   *  by their greatest common divisor.
   */
  public RationalNumber canonicalize() {
    BigInteger gcd = numerator.gcd(denominator);
    if (! BigInteger.ONE.equals(gcd)) {
      numerator = numerator.divide(gcd);
      denominator =  denominator.divide(gcd);
    }
    return this;
  }

  private void adjustSigns() {
    if (denominator.signum() == -1) {
      // Always keep denominator positive
      denominator = denominator.negate();
      numerator = numerator.negate();
    }
  }

  public static void main(String[] args) throws Exception {
    RationalNumber r1 = new RationalNumber(BigInteger.valueOf(20),
                                           BigInteger.valueOf(-5));
    System.out.println("r1 = " + r1);
    RationalNumber r2 = new RationalNumber(BigInteger.valueOf(2),
                                           BigInteger.valueOf(-5));
    System.out.println("r2 = " + r2);
    System.out.println("r1 + r2 = " + r1.add(r2).canonicalize());
    System.out.println("r1 - r2 = " + r1.subtract(r2).canonicalize());
    System.out.println("r2 - r1 = " + r2.subtract(r1).canonicalize());
    System.out.println("r1*r2 = " + r1.multiply(r2).canonicalize());
    System.out.println("r1/r2 = " + r1.divide(r2).canonicalize());
    System.out.println("r2/r1 = " + r2.divide(r1).canonicalize());
  }

}
