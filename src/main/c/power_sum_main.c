#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bernoulli_power_sum.h"
#include "central_factorial_power_sum.h"
#include "euler_power_sum.h"
#include "faulhaber_power_sum.h"
#include "power_sum.h"
#include "stirling_power_sum.h"

static char * VALID_OPTIONS[] = {"-c", "-f", "-h", "-s", "-sv"};

static void usage(char * commandName) {
  fprintf(stderr, "Usage: %s"
  " (-c|-f|-h|-s|-sv) [<power>] [<numTerms>]\n\n"
  "<power> and <numTerms> should be greater than or equal to 0\n\n"
  "Examples:\n"
  "To print the help on usage:\n"
  "%s -h or just %s\n\n"
  "To print coefficients in the formula for power 10:\n"
  "%s -c 10\n\n"
  "To print formula for sum for power 5:\n"
  "%s -f 5\n\n"
  "To print the sum of series for power 6 for the first 20 terms:\n"
  "The sum will be computed using the formula\n"
  "%s -s 6 20\n\n"
  "To compute the sum in two ways one using the formula and the\n"
  "other with actual series expansion "
  "and verify the results for correctness\n"
  "%s -sv 6 20\n\n"
  "If <numTerms> is missing, a default of 20 is assumed\n", commandName,
  commandName, commandName, commandName, commandName, commandName, commandName);
}

static void error(char * commandName, char * message) {
  fprintf(stderr, "%s\n", message);
  usage(commandName);
  exit(EXIT_FAILURE);
}

static void error2(char * commandName, char * message1, char * message2) {
  fprintf(stderr, "%s%s\n", message1, message2);
  usage(commandName);
  exit(EXIT_FAILURE);
}

static unsigned long isValidOption(char * option) {
  unsigned long i;
  for (i = 0; i < sizeof(VALID_OPTIONS)/sizeof(char *); i++) {
    if (strcmp(VALID_OPTIONS[i], option) == 0) {
      return 1UL;
    }
  }
  return 0UL;
}

static void printCpuTime(struct timespec * before, struct timespec * after) {
  printf("Time taken = ");
  long cpuNano = (after->tv_sec - before->tv_sec)*1000000000L
                      + after->tv_nsec - before->tv_nsec;
  printf("%ld\n", cpuNano);
}

static void printFaulhaberTitle() {
  printf("Faulhaber: -------------------------------------\n");
}

static void printBernoulliTitle() {
  printf("Bernoulli: -------------------------------------\n");
}

static void printStirlingTitle() {
  printf("Stirling: --------------------------------------\n");
}

static void printEulerTitle() {
  printf("Euler: -----------------------------------------\n");
}

static void printCentralFactorialTitle() {
  printf("Central Factorial: -----------------------------\n");
}

static void printCoefficients(mpq_t * coeffs, long numCoeffs) {
  long t;

  for (t = 0; t < numCoeffs; t++) {
    gmp_printf(" %Qd", coeffs[t]);
  }
  printf("\n");
}

static void getCoefficientsTimed(PowerSum * psPtr, long power) {
  struct timespec before;
  struct timespec after;
  mpq_t * coeffs;
  long numCoeffs;

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
  psPtr->getCoefficients(power, &coeffs, &numCoeffs);
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
  printCoefficients(coeffs, numCoeffs);
  printCpuTime(&before, &after);
  psPtr->freeCoefficients(coeffs, numCoeffs);
}

static void computeAndPrintSumTimed(PowerSum * psPtr, long power, long numTerms,
                                    mpz_t sum) {
  long coeffInitTime;
  long sumTime;

  psPtr->computeSumWithTimeStat(power, numTerms, sum, &coeffInitTime, &sumTime);
  gmp_printf("Sum computed = %Zd\n", sum);
  printf("Time taken = ");
  long total = coeffInitTime + sumTime;
  printf("%ld:%ld:%ld\n", total, coeffInitTime, sumTime);
}

int main(int argc, char ** argv) {

  /* Validate command arguments */
  if (argc <= 1 || strcmp(argv[1], "-h") == 0) {
    usage(argv[0]);
    return EXIT_SUCCESS;
  }

  if (isValidOption(argv[1]) == 0) {
    error2(argv[0], "Invalid option: ", argv[1]);
  }

  if (argc < 3) {
    error(argv[0], "Missing mandatory argument: <power>");
  }

  long power = 1;
  long numTerms = 20;
  power = atol(argv[2]);
  if (power < 0) {
    error2(argv[0], "Invalid power: ", argv[2]);
  }
  if (argc > 3) {
    numTerms = atol(argv[3]);
    if (numTerms < 0) {
      error2(argv[0], "Invalid number of terms: ", argv[3]);
    }
  }

  /* Initialize function pointers */
  PowerSum fps;
  initFaulhaberPowerSum(&fps);
  PowerSum bps;
  initBernoulliPowerSum(&bps);
  PowerSum sps;
  initStirlingPowerSum(&sps);
  PowerSum eps;
  initEulerPowerSum(&eps);
  PowerSum cps;
  initCentralFactorialPowerSum(&cps);
  struct timespec before;
  struct timespec after;

  /* Do what the user asked */
  if (strcmp(argv[1], "-c") == 0) {
    printf("Computing coefficients for power %ld\n", power);
    printFaulhaberTitle();
    getCoefficientsTimed(&fps, power);
    printBernoulliTitle();
    getCoefficientsTimed(&bps, power);
    printStirlingTitle();
    getCoefficientsTimed(&sps, power);
    printEulerTitle();
    getCoefficientsTimed(&eps, power);
    printCentralFactorialTitle();
    getCoefficientsTimed(&cps, power);
  } else if (strcmp(argv[1], "-f") == 0) {
    printFaulhaberTitle();
    fps.printSumFormula(power, stdout);
    printBernoulliTitle();
    bps.printSumFormula(power, stdout);
    printStirlingTitle();
    sps.printSumFormula(power, stdout);
    printEulerTitle();
    eps.printSumFormula(power, stdout);
    printCentralFactorialTitle();
    cps.printSumFormula(power, stdout);
  } else {
    mpz_t sumFaulhaber;
    mpz_t sumBernoulli;
    mpz_t sumStirling;
    mpz_t sumEuler;
    mpz_t sumCentral;

    mpz_init(sumFaulhaber);
    mpz_init(sumBernoulli);
    mpz_init(sumStirling);
    mpz_init(sumEuler);
    mpz_init(sumCentral);
    printf("Computing S(%ld, %ld)\n", power, numTerms);
    printFaulhaberTitle();
    computeAndPrintSumTimed(&fps, power, numTerms, sumFaulhaber);
    printBernoulliTitle();
    computeAndPrintSumTimed(&bps, power, numTerms, sumBernoulli);
    printStirlingTitle();
    computeAndPrintSumTimed(&sps, power, numTerms, sumStirling);
    printEulerTitle();
    computeAndPrintSumTimed(&eps, power, numTerms, sumEuler);
    printCentralFactorialTitle();
    computeAndPrintSumTimed(&cps, power, numTerms, sumCentral);
    if (strcmp(argv[1], "-sv") == 0) {
      mpz_t sumSeries;
      mpz_init(sumSeries);
      printf("Series addition:--------------------------------\n");
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &before);
      sps.computeSumUsingSeries(power, numTerms, sumSeries);
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &after);
      gmp_printf("Sum computed = %Zd\n", sumSeries);
      printCpuTime(&before, &after);
      if (mpz_cmp(sumFaulhaber, sumSeries) == 0){
        printf("The sum matches with Faulhaber formula :-)\n");
      } else {
        printf("The sums do not match for Faulhaber formula :-(\n");
      }
      if (mpz_cmp(sumBernoulli, sumSeries) == 0){
        printf("The sum matches with Bernoulli formula :-)\n");
      } else {
        printf("The sums do not match for Bernoulli formula :-(\n");
      }
      if (mpz_cmp(sumStirling, sumSeries) == 0){
        printf("The sum matches with Stirling formula :-)\n");
      } else {
        printf("The sums do not match for Stirling formula :-(\n");
      }
      if (mpz_cmp(sumEuler, sumSeries) == 0){
        printf("The sum matches with Euler formula :-)\n");
      } else {
        printf("The sums do not match for Euler formula :-(\n");
      }
      if (mpz_cmp(sumCentral, sumSeries) == 0){
        printf("The sum matches with Central Factorial formula :-)\n");
      } else {
        printf("The sums do not match for Central Factorial formula :-(\n");
      }
      mpz_clear(sumSeries);
    }
    mpz_clear(sumCentral);
    mpz_clear(sumEuler);
    mpz_clear(sumStirling);
    mpz_clear(sumBernoulli);
    mpz_clear(sumFaulhaber);
  }
  return EXIT_SUCCESS;
}
