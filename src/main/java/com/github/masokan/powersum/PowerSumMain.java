package com.github.masokan.powersum;

import java.math.BigInteger;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class PowerSumMain {
  private static final String USAGE = 
    "Usage: %command% (-c|-f|-s|-sv) [<power>] [<numTerms>]\n\n" +
    "<power> and <numTerms> should be greater than or equal to 0\n\n" +
    "Examples:\n" +
    "To print the help on usage:\n" +
    "%command% -h or just %command%\n\n" +
    "To print coefficients in the formula for power 10:\n" +
    "%command% -c 10\n\n" +
    "To print formula for sum for power 5:\n" +
    "%command% -f 5\n\n" +
    "To print the sum of series for power 6 for the first 20 terms:\n" +
    "The sum will be computed using the formula\n" +
    "%command% -s 6 20\n\n" +
    "To compute the sum in two ways one using the formula and the\n" +
    "other with actual series expansion " +
    "and verify the results for correctness\n" +
    "%command% -sv 6 20\n\n" +
    "If <numTerms> is missing, a default of 20 is assumed";

  private static final String[] VALID_OPTIONS
      = {"-c", "-f", "-h", "-s", "-sv"};

  private static void usage(String commandName) {
    System.err.println(USAGE.replaceAll("%command%", commandName));
  }

  private static void error(String commandName, String message) {
    System.err.println(message);
    usage(commandName);
    System.exit(1);
  }

  private static void printCoefficients(List<RationalNumber> coeffs) {
    for (int t = 0; t < coeffs.size(); t++) {
      System.out.print(" " + coeffs.get(t));
    }
    System.out.println();
  }

  public static void printCpuTime(long before, long after) {
    System.out.print("Time taken = ");
    long elapsedNano = after - before;
    System.out.println(elapsedNano);
  }

  private static void printFaulhaberTitle() {
    System.out.println("Faulhaber: -------------------------------------");
  }

  private static void printBernoulliTitle() {
    System.out.println("Bernoulli: -------------------------------------");
  }

  private static void printStirlingTitle() {
    System.out.println("Stirling: --------------------------------------");
  }

  private static void printEulerTitle() {
    System.out.println("Euler: -----------------------------------------");
  }

  private static void printCentralFactorialTitle() {
    System.out.println("Central Factorial: -----------------------------");
  }

  private static void getCoefficientsTimed(PowerSum ps, long power) {
    long before = PowerSum.getThreadCpuTime();
    List<RationalNumber> coeffs = ps.getCoefficients(power);
    long after = PowerSum.getThreadCpuTime();
    printCoefficients(coeffs);
    printCpuTime(before, after);
  }
  
  private static BigInteger computeAndPrintSumTimed(PowerSum ps, long power,
                                                    long numTerms) {
    long[] stat = new long[2];
    BigInteger sum = ps.computeSumWithTimeStat(power, numTerms, stat);
    long totalTime = stat[0] + stat[1];
    System.out.println("Sum computed = " + sum);
    System.out.print("Time taken = ");
    System.out.println(totalTime + ":" + stat[0] + ":" + stat[1]);
    return sum;
  }

  public static void main(String[] args) {
    HashSet<String> validOptions = new HashSet<>(Arrays.asList(VALID_OPTIONS));

    PowerSum fps = new FaulhaberPowerSum();
    PowerSum bps = new BernoulliPowerSum();
    PowerSum sps = new StirlingPowerSum();
    PowerSum eps = new EulerPowerSum();
    PowerSum cps = new CentralFactorialPowerSum();
    String commandName = "java -jar PowerSum.jar";
    // Validate command arguments
    if (args.length == 0 || "-h".equals(args[0])) {
      usage(commandName);
      System.exit(0);
    }

    if (! validOptions.contains(args[0])) {
      error(commandName, "Invalid option: " + args[0]);
    }

    if (args.length < 2) {
      error(commandName, "Missing mandatory argument: <power>");
    }

    long power = 1;
    long numTerms = 20;
    try {
      power = Long.parseLong(args[1]);
      if (power < 0) {
        throw new NumberFormatException("");
      }
    } catch(NumberFormatException nfe) {
      error(commandName, "Invalid power: " + args[1]);
    }
    if (args.length > 2) {
      try {
        numTerms = Long.parseLong(args[2]);
        if (numTerms < 0) {
          throw new NumberFormatException("");
        }
      } catch(NumberFormatException nfe) {
        error(commandName, "Invalid number of terms: " + args[2]);
      }
    }

    long before;
    long after;
    // Do what the user asked
    if ("-c".equals(args[0])) {
      System.out.println("Computing coefficients for power " + power);
      printFaulhaberTitle();
      getCoefficientsTimed(fps, power);
      printBernoulliTitle();
      getCoefficientsTimed(bps, power);
      printStirlingTitle();
      getCoefficientsTimed(sps, power);
      printEulerTitle();
      getCoefficientsTimed(eps, power);
      printCentralFactorialTitle();
      getCoefficientsTimed(cps, power);
    } else if ("-f".equals(args[0])) {
      printFaulhaberTitle();
      fps.printSumFormula(power, System.out);
      printBernoulliTitle();
      bps.printSumFormula(power, System.out);
      printStirlingTitle();
      sps.printSumFormula(power, System.out);
      printEulerTitle();
      eps.printSumFormula(power, System.out);
      printCentralFactorialTitle();
      cps.printSumFormula(power, System.out);
    } else {
      System.out.println("Computing S(" + power + ", " + numTerms + ")");
      printFaulhaberTitle();
      BigInteger sumFaulhaber = computeAndPrintSumTimed(fps, power, numTerms);
      printBernoulliTitle();
      BigInteger sumBernoulli = computeAndPrintSumTimed(bps, power, numTerms);
      printStirlingTitle();
      BigInteger sumStirling = computeAndPrintSumTimed(sps, power, numTerms);
      printEulerTitle();
      BigInteger sumEuler = computeAndPrintSumTimed(eps, power, numTerms);
      printCentralFactorialTitle();
      BigInteger sumCentral = computeAndPrintSumTimed(cps, power, numTerms);
      if ("-sv".equals(args[0])) {
        System.out.println("Series addition:--------------------------------");
        before = PowerSum.getThreadCpuTime();
        BigInteger sumFromSeries = sps.computeSumUsingSeries(power, numTerms);
        after = PowerSum.getThreadCpuTime();
        System.out.println("Sum computed = " + sumFromSeries);
        printCpuTime(before, after);
        if (sumFaulhaber.equals(sumFromSeries)) {
          System.out.println("The sum matches with Faulhaber formula :-)");
        } else {
          System.out.println("Sums do not match for Faulhaber formula :-(");
        }
        if (sumBernoulli.equals(sumFromSeries)) {
          System.out.println("The sum matches with Bernoulli formula :-)");
        } else {
          System.out.println("Sums do not match for Bernoulli formula :-(");
        }
        if (sumStirling.equals(sumFromSeries)) {
          System.out.println("The sum matches with Stirling formula :-)");
        } else {
          System.out.println("Sums do not match for Stirling formula :-(");
        }
        if (sumEuler.equals(sumFromSeries)) {
          System.out.println("The sum matches with Euler formula :-)");
        } else {
          System.out.println("Sums do not match for Euler formula :-(");
        }
        if (sumCentral.equals(sumFromSeries)) {
          System.out.println("The sum matches with Central Factorial formula :-)");
        } else {
          System.out.println("Sums do not match for Central Factorial formula :-(");
        }
      }
    }
    System.exit(0);
  }

}
