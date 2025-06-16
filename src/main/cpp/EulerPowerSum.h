#ifndef EULER_POWERSUM_H
#define EULER_POWERSUM_H

#include "PowerSum.h"

class EulerPowerSum : public PowerSum {
  public:
    EulerPowerSum();
    virtual vector<mpq_class> getCoefficients(long power);
    virtual void printSumFormula(long power, ostream &out);
    virtual mpz_class computeSumWithTimeStat(long power, long n,
                                             vector<long> & stat);
    virtual ~EulerPowerSum();

  private:
    vector<mpz_class> getCoefficients(long power, long maxNumCoefficients);
    void printTerm(long start, long numTerms, ostream & out);

};

#endif
