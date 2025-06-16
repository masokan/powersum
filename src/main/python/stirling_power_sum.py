from gmpy2 import mpz, mpq
from power_sum import PowerSum
import time

class StirlingPowerSum(PowerSum):

  def __init__(self) -> None:
    return

  def get_coefficients(self, power: int) -> list[mpq]:
    coeffs: list[mpz] = self.__get_coefficients(power, power + 1)
    return self.to_rational(coeffs)

  def print_sum_formula(self, power: int) -> None:
    if power < 0:
      return

    coeffs: list[mpz] = self.__get_coefficients(power, power + 1)
    first_time: bool = True
    for t in range(len(coeffs)):
      coeff: mpz = coeffs[t];
      if (coeff != 0):
        if first_time:
          print('   ', end='')
          first_time = False
        else:
          print(' + ', end='')
        if coeff > 1:
          # Coefficient is greater than 1
          print(f'{coeff}', end='')
        self.__print_factors(t)
        if (t > 0):
          print(f'/{(t + 1)}', end='')
    print()

  def compute_sum_with_time_stat(self, power: int, n: int,
                                 stat: list[int]) -> mpz:
    sum: mpz = 0
    stat.clear()

    if (power < 0 or n < 0):
      stat.append(0)
      stat.append(0)
      return sum

    before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    coeffs: list[mpz] = self.__get_coefficients(power, n + 1)
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)

    before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    falling_factorial: mpz = n + 1
    factor: mpz = 0
    num_terms_to_compute = len(coeffs)
    if num_terms_to_compute > (n + 1):
      num_terms_to_compute = n + 1
    for t in range(num_terms_to_compute):
      # We know that following division is exact and we do this first before
      # multiplication to avoid generating a large intermediate value
      factor = falling_factorial // (t + 1)
      sum = sum + coeffs[t]*factor
      falling_factorial = falling_factorial * (n - t)
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)
    return sum

   # The coefficients are Stirling numbers of second kind which are defined as
   # below:
   # S(0, 0) = 1
   # S(m, 0) = 0 for all m > 0
   # S(m, j) = 0 for all j > m
   # S(m, j) = S(m-1, j-1) + j*S(m-1, j)
   # m is the power of the series (m >= 0)
   # j is the index of the term in the sum formula (j >= 0)
  def __get_coefficients(self, power: int, max_coefficients: int) -> list[mpz]:
    coeffs: list[mpz] = []

    if power < 0:
      return coeffs

    if power > 0:
      Sm_1j_1: mpz = 1 # S(0, 0)
      Smj: mpz = 0

      # Initialize all coefficients to 0
      for i in range(power + 1):
        coeffs.append(0)

      for current_power in range(1, power + 1):
        term_range = min(power + 1, max_coefficients + 1)
        # Initialize other coefficients for the current power
        for term in range(1, term_range):
          Smj = Sm_1j_1 + coeffs[term]*term
          Sm_1j_1 = coeffs[term]
          coeffs[term] = Smj
        Sm_1j_1 = coeffs[0]
    else:
      coeffs.append(1)

    return coeffs

  def __print_factors(self, term: int) -> None:
    for i in range(term + 1):
      match i:
        case 0:
          print('(n + 1)', end='')
        case 1:
          print('n', end='')
        case _:
          print(f'(n - {(i - 1)})', end='')

