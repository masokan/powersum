#ifndef BERNOULLI_POWERSUM_H
#define BERNOULLI_POWERSUM_H

#include "PowerSum.h"

class BernoulliPowerSum : public PowerSum {
  public:
    BernoulliPowerSum();
    virtual vector<mpq_class> getCoefficients(long power);
    virtual void printSumFormula(long power, ostream &out);
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat);
    virtual ~BernoulliPowerSum();

  private:
    mpq_class computeNextCoefficient(vector<mpq_class> & current, long m);

};

#endif
