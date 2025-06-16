#include <stdlib.h>
#include <time.h>

#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "BernoulliPowerSum.h"
#include "CentralFactorialPowerSum.h"
#include "EulerPowerSum.h"
#include "FaulhaberPowerSum.h"
#include "PowerSum.h"
#include "StirlingPowerSum.h"

using std::cerr;
using std::cout;
using std::endl;
using std::invalid_argument;
using std::set;
using std::stol;
using std::string;

static void usage(string & commandName) {
  cerr << "Usage: " << commandName
  << " (-c|-f|-h|-s|-sv) [<power>] [<numTerms>]" << endl << endl
  << "<power> and <numTerms> should be greater than or equal to 0" << endl
  << endl << "Examples:\n"
  << "To print the help on usage:" << endl
  << commandName << " -h or just " << commandName << endl << endl
  << "To print coefficients in the formula for power 10:" << endl
  << commandName << " -c 10" << endl << endl
  << "To print formula for sum for power 5:" << endl
  << commandName << " -f 5" << endl << endl
  << "To print the sum of series for power 6 for the first 20 terms:" << endl
  << "The sum will be computed using the formula" << endl
  << commandName << " -s 6 20" << endl << endl
  << "To compute the sum in two ways one using the formula and the" << endl
  << "other with actual series expansion "
  << "and verify the results for correctness" << endl
  << commandName << " -sv 6 20" << endl << endl
  << "If <numTerms> is missing, a default of 20 is assumed" << endl;
}

static void error(string & commandName, string message) {
  cerr << message << endl;
  usage(commandName);
  exit(EXIT_FAILURE);
}

static void printCpuTime(struct timespec & before, struct timespec & after) {
  cout << "Time taken = ";
  long cpuNano = (after.tv_sec - before.tv_sec)*1000000000L
                      + after.tv_nsec - before.tv_nsec;
  cout << cpuNano << endl;
}

static void printCoefficients(vector<mpq_class> & coeffs, ostream &out) {
  for (size_t t = 0; t < coeffs.size(); t++) {
    out << " " << coeffs[t];
  }
  out << endl;
}

static void printFaulhaberTitle() {
  cout << "Faulhaber: -------------------------------------" << endl;
}

static void printBernoulliTitle() {
  cout << "Bernoulli: -------------------------------------" << endl;
}

static void printStirlingTitle() {
  cout << "Stirling: --------------------------------------" << endl;
}

static void printEulerTitle() {
  cout << "Euler: -----------------------------------------" << endl;
}

static void printCentralFactorialTitle() {
  cout << "Central Factorial: -----------------------------" << endl;
}

static void getCoefficientsTimed(PowerSum & ps, long power) {
  struct timespec before;
  struct timespec after;

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  vector<mpq_class> coeffs = ps.getCoefficients(power);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  printCoefficients(coeffs, cout);
  printCpuTime(before, after);
}

static mpz_class computeAndPrintSumTimed(PowerSum & ps, long power,
                                         long numTerms) {
  vector<long> stat;

  mpz_class sum = ps.computeSumWithTimeStat(power, numTerms, stat);
  cout << "Sum computed = " << sum << endl;
  cout << "Time taken = ";
  long total = stat[0] + stat[1];
  cout << total << ':' << stat[0] << ':' << stat[1] << endl;
  return sum;
}

int main(int argc, char ** argv) {
  set<string> validOptions;
  validOptions.insert("-c");
  validOptions.insert("-f");
  validOptions.insert("-h");
  validOptions.insert("-s");
  validOptions.insert("-sv");
  struct timespec before;
  struct timespec after;

  string * args = new string[argc];
  for (int i = 0; i < argc; i++) {
    args[i] = argv[i];
  }

  // Validate command arguments
  if (argc <= 1 || args[1] == "-h") {
    usage(args[0]);
    return EXIT_SUCCESS;
  }

  if (validOptions.find(args[1]) == validOptions.end()) {
    error(args[0], string("Invalid option: " + args[1]));
  }

  if (argc < 3) {
    error(args[0], "Missing mandatory argument: <power>");
  }

  long power = 1;
  long numTerms = 20;
  try {
    power = stol(args[2]);
    if (power < 0) {
      throw invalid_argument("");
    }
  } catch(...) {
    error(args[0], "Invalid power: " + args[2]);
  }
  if (argc > 3) {
    try {
      numTerms = stol(args[3]);
      if (numTerms < 0) {
        throw invalid_argument("");
      }
    } catch(...) {
      error(args[0], "Invalid number of terms: " + args[3]);
    }
  }

  // Do what the user asked
  FaulhaberPowerSum fps;
  BernoulliPowerSum bps;
  StirlingPowerSum sps;
  EulerPowerSum eps;
  CentralFactorialPowerSum cps;
  mpz_class sumFaulhaber;
  mpz_class sumBernoulli;
  mpz_class sumStirling;
  mpz_class sumEuler;
  mpz_class sumCentral;

  if (args[1] == "-c") {
    cout << "Computing coefficients for power " << power << endl;
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
  } else if (args[1] == "-f") {
    printFaulhaberTitle();
    fps.printSumFormula(power, cout);
    printBernoulliTitle();
    bps.printSumFormula(power, cout);
    printStirlingTitle();
    sps.printSumFormula(power, cout);
    printEulerTitle();
    eps.printSumFormula(power, cout);
    printCentralFactorialTitle();
    cps.printSumFormula(power, cout);
  } else {
    cout << "Computing S(" << power << ", " << numTerms << ")" << endl;
    printFaulhaberTitle();
    sumFaulhaber = computeAndPrintSumTimed(fps, power, numTerms);
    printBernoulliTitle();
    sumBernoulli = computeAndPrintSumTimed(bps, power, numTerms);
    printStirlingTitle();
    sumStirling = computeAndPrintSumTimed(sps, power, numTerms);
    printEulerTitle();
    sumEuler = computeAndPrintSumTimed(eps, power, numTerms);
    printCentralFactorialTitle();
    sumCentral = computeAndPrintSumTimed(cps, power, numTerms);
    if (args[1] == "-sv") {
      cout << "Series addition:--------------------------------" << endl;
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
      mpz_class sumFromSeries = sps.computeSumUsingSeries(power, numTerms);
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
      cout << "Sum computed = " << sumFromSeries << endl;
      printCpuTime(before, after);
      if (sumFaulhaber == sumFromSeries) {
        cout << "The sum matches with Faulhaber formula :-)" << endl;
      } else {
        cout << "The sums do not match for Faulhaber formula :-(" << endl;
      }
      if (sumBernoulli == sumFromSeries) {
        cout << "The sum matches with Bernoulli formula :-)" << endl;
      } else {
        cout << "The sums do not match for Bernoulli formula :-(" << endl;
      }
      if (sumStirling == sumFromSeries) {
        cout << "The sum matches with Stirling formula :-)" << endl;
      } else {
        cout << "The sums do not match for Stirling formula :-(" << endl;
      }
      if (sumEuler == sumFromSeries) {
        cout << "The sum matches with Euler formula :-)" << endl;
      } else {
        cout << "The sums do not match for Euler formula :-(\n";
      }
      if (sumCentral == sumFromSeries) {
        cout << "The sum matches with Central Factorial formula :-)" << endl;
      } else {
        cout << "The sums do not match for Central Factorial formula :-("
             << endl;
      }
    }
  }
  delete[] args;
  return EXIT_SUCCESS;
}
