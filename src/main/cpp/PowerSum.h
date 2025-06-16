#ifndef POWERSUM_H
#define POWERSUM_H

#include <time.h>
#include <gmpxx.h>

#include <ostream>
#include <vector>

using std::ostream;
using std::vector;

class PowerSum {
  public:
    PowerSum() {}
    /* To get the coefficients in the sum formula.
     * Parameters:
     *   power - desired power (IN)
     * Return value
     *   a vector of coefficients
     */
    virtual vector<mpq_class> getCoefficients(long power) = 0;

    /* To print the sum formula for a given power
     * Parameters:
     *   power - desired power (IN)
     *   out - output stream where the formula should be printed (IN)
     */
    virtual void printSumFormula(long power, ostream &out) = 0;

    /* To compute the sum for a specific power and the number of terms
     * Parameters:
     *   power - desired power (IN)
     *   numTerms - number of terms (IN)
     * Return value
     *   sum computed
     */
    virtual mpz_class computeSum(long power, long n);

    /* To compute the sum for a specific power and the number of terms and
     * obtain CPU times
     * Parameters:
     *   power - desired power (IN)
     *   numTerms - number of terms (IN)
     *   stat - vector in which the coefficient initialization time and
     *          summation time (in nanosec.) are returned (IN/OUT)
     * Return value
     *   sum computed
     */
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat) = 0;
    /* To compute the sum for a specific power and the number of terms using the
     * simple implementation (series summation)
     * Parameters:
     *   power - desired power (IN)
     *   numTerms - number of terms (IN)
     * Return value
     *   sum computed
     */
    virtual mpz_class computeSumUsingSeries(long power, long n);

    virtual ~PowerSum() {}
    // Some useful implementations for use in derived classes
  protected:
    long computeCpuTime(struct timespec & before, struct timespec & after);
    mpz_class nCr(long n, long r);
};

#endif
