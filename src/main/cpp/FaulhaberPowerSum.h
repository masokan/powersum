#ifndef FAULHABER_POWERSUM_H
#define FAULHABER_POWERSUM_H

#include "PowerSum.h"

class FaulhaberPowerSum : public PowerSum {
  public:
    FaulhaberPowerSum();
    virtual vector<mpq_class> getCoefficients(long power);
    virtual void printSumFormula(long power, ostream &out);
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat);
    virtual ~FaulhaberPowerSum();

  private:
    mpz_class postProcessRow(vector<mpz_class> & row);
    mpz_class createRowForEvenPower(long nLimit, long rowNum,
                                    vector<mpz_class> & row);
    mpz_class createRowForOddPower(long nLimit, long rowNum,
                                   vector<mpz_class> & row);
};

#endif
