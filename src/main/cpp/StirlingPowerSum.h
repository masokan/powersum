#ifndef STIRLING_POWERSUM_H
#define STIRLING_POWERSUM_H

#include "PowerSum.h"

class StirlingPowerSum : public PowerSum {
  public:
    StirlingPowerSum();
    virtual vector<mpq_class> getCoefficients(long power);
    virtual void printSumFormula(long power, ostream &out);
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat);
    virtual ~StirlingPowerSum();

  private:
    vector<mpz_class> getCoefficients(long power, long maxNumCoefficients);
    void printFactors(long term, ostream & out);

};

#endif
