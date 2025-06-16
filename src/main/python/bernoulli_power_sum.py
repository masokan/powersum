from gmpy2 import mpz, mpq
from power_sum import PowerSum
import time

class BernoulliPowerSum(PowerSum):

  def __init__(self) -> None:
    return

  # Bernoulli numbers are defined by the following recurrence relation:
  # B(0) = 1
  # B(m) = -(Binom((m + 1), 0)B(0) + Binom((m + 1), 1)B(1)
  #          + ... + Binom((m + 1), (m - 1))B(m - 1))
  # where Binom(i, j) is the binomial coefficient which evaluates to:
  # i!/{(i - j)!j!}
  def get_coefficients(self, power: int) -> list[mpq]:
    coeffs: list[mpq] = []
    if power < 0:
      return coeffs

    i: int

    for i in range(0, power + 1):
     coeffs.append(0)

    # Initialize B(0)
    coeffs[0] = mpq(1)
    if power > 0:
      # Initialize B(1)
      coeffs[1] = mpq(-1, 2)

    for i in range(2, power + 1):
      if ((i & 1) == 1):
        # Odd coefficients above 1 are 0s
        coeffs[i] = mpq(0)
      else:
        coeffs[i] = self.__compute_next_coefficient(coeffs, i)

    return coeffs

  def print_sum_formula(self, power: int) -> None:
    if power < 0:
      return

    coeffs: list[mpq] = self.get_coefficients(power)
    binom: mpz = 1
    binomN: int = power + 1
    binomR: int = 1
    p: int = power + 1

    print('{ ', end='')
    for i in range(0, len(coeffs)):
      if ((i & 1) == 0 or i == 1):
        if (i != 0):
          print(' + ', end='')
        if (coeffs[i] != 1):
          print(f'({coeffs[i]})', end='')
        if (binom != 1):
          print(binom, end='')
        print('(n + 1)', end='')
        if (p != 1):
          print(f'^{p}', end='')
      p -= 1
      binom = binom*binomN//binomR
      binomN -= 1
      binomR += 1
    
    print(' }', end='')
    if (power > 0):
      print(f'/{(power + 1)}', end='')
    print()

  def compute_sum_with_time_stat(self, power: int, n: int,
                                 stat: list[int]) -> mpz:
    stat.clear()

    if (power < 0 or n < 0):
      stat.append(0)
      stat.append(0)
      result: mpz = 0
      return result

    before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    coeffs: list[mpq] = self.get_coefficients(power)
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)

    before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    sum: mpq = 0
    p: mpz = n + 1
    binomN: int = power + 1
    binomR: int = 1
    binom: mpz = binomN

    # We compute the terms in reverse order to avoid an additional division
    # operation in the loop
    for i in range(len(coeffs) - 1, -1, -1):
      if ((i & 1) == 0 or i == 1):
        sum += binom*coeffs[i]*p
      p *= (n + 1)
      binomN -= 1
      binomR += 1
      binom = binom*binomN//binomR

    sum /= (power + 1)
    (numerator, denominator) = sum.as_integer_ratio()
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)
    return numerator

  def __compute_next_coefficient(self, current: list[mpq], m: int) -> mpq:
    binom: mpz = 1
    coeff: mpq = 0
    k: int

    for k in range (0, m):
      if (((k & 1) == 0) or k == 1):
        # Even coefficient or the first odd one
        coeff += current[k]*binom
      binom *= (m + 1 - k)
      binom //= (k + 1)
    coeff = -coeff/binom
    return coeff
