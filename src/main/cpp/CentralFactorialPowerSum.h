#ifndef CENTRAL_FACTORIAL_POWER_SUM_H
#define CENTRAL_FACTORIAL_POWER_SUM_H

#include "PowerSum.h"

class CentralFactorialPowerSum : public PowerSum {
  public:
    CentralFactorialPowerSum();
    virtual vector<mpq_class> getCoefficients(long power);
    virtual void printSumFormula(long power, ostream &out);
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat);
    virtual ~CentralFactorialPowerSum();

  private:
    vector<mpz_class> getCoefficients(long power, long maxN);
    void printFallingFactorial(long start, long numTerms, ostream & out);

};

#endif
