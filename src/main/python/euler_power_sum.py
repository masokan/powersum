from gmpy2 import mpz, mpq
from power_sum import PowerSum
import time

class EulerPowerSum(PowerSum):

  def __init__(self) -> None:
    return

  def get_coefficients(self, power: int) -> list[mpq]:
    coeffs: list[mpz] = self.__get_coefficients(power, power)
    return self.to_rational(coeffs)

  def print_sum_formula(self, power: int) -> None:
    if power < 0:
      return

    coeffs: list[mpz] = self.__get_coefficients(power, power)

    # Compute (power + 1)!
    power_factorial: mpz = 1
    for i in range(2, power + 2):
      power_factorial = power_factorial*i

    print('{ ', end='')
    for j in range(0, len(coeffs)):
      if (coeffs[j] > 1):
        # Coefficient is greater than 1
        print(coeffs[j], end='')
      if (coeffs[j] > 0):
        self.__print_term(j, power + 1)
      if (j < (len(coeffs) - 2)):
        print(' + ', end='')
    print(' }', end='')
    if (power_factorial != 1):
      print(f'/{power_factorial}', end='')
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
    coeffs: list[mpz] = self.__get_coefficients(power, n)
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)

    before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    temp: mpz
    j: int
    first_time: bool = True

    for j in range(0, len(coeffs)):
      if (n + j >= power):
        if (first_time):
          first_time = False
          temp = self.nCr(n + j + 1, power + 1)
        else:
          temp *= (n + j + 1)
          temp //= (n + j - power)
        sum += coeffs[j]*temp

    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)
    return sum

   # The Euler numbers of first kind are as defined below:
   # E(i, j) = 1 when j = 0
   # E(i, j) = (j + 1)E(n - 1, j) + (n - j)E(n - 1, j - 1)
  def __get_coefficients(self, power: int, max_n: int) -> list[mpz]:
    coeffs: list[mpz] = []

    if power < 0:
      return coeffs

    temp: mpz
    E_i_1_j_1: mpz
    i: int
    j: int
    k: int

    coeffs.append(1)
    for k in range(1, power + 1):
      coeffs.append(0)
    for i in range(1, power + 1):
      odd_power = ((i & 1) == 1)
      halfI = i >> 1
      if odd_power:
        halfLimit = halfI
      else:
        halfLimit = (halfI - 1)
      limit = halfLimit
      # No need to initialize more than maximum n since falling factorial in
      # other terms will be 0
      if (halfLimit > max_n):
        limit = max_n
      for j in range(0, limit + 1):
        if (j == 0):
          E_i_1_j_1 = coeffs[j]
          coeffs[j] = 1
        else:
          temp = (j + 1)*coeffs[j] + (i - j)*E_i_1_j_1
          E_i_1_j_1 = coeffs[j]
          coeffs[j] = temp

      j = halfLimit
      # Since coefficients are mirror reflection w.r.t. the central point,
      # we can reverse copy whatever we initialized so far
      if (odd_power):
        k = 1
        while(j + k < i):
          coeffs[j + k] = coeffs[j - k]
          k += 1
      else:
        k = 1
        while(j + k < i):
          coeffs[j + k] = coeffs[j - k + 1]
          k += 1
      coeffs[i] = 0
    return coeffs

  def __print_term(self, start: int, num_terms: int) -> None:
    start += 1
    for j in range(0, num_terms):
      if (start > 0):
        print(f'(n + {start})', end='')
      elif (start == 0):
        print('n', end='')
      else:
        print(f'(n - {-start})', end='')
      start -= 1
