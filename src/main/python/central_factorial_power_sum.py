from gmpy2 import mpz, mpq
from power_sum import PowerSum
import time

class CentralFactorialPowerSum(PowerSum):

  def __init__(self) -> None:
    return

  def get_coefficients(self, power: int) -> list[mpq]:
    coeffs: list[mpz] = self.__get_coefficients(power, power)
    return self.to_rational(coeffs)

  def print_sum_formula(self, power: int) -> None:
    if power < 0:
      return

    even_power: bool = ((power & 1) == 0)

    if power > 0:
      coeffs: list[mpz] = self.__get_coefficients(power, power)
      numCoeffs = len(coeffs)
      # Coefficient 0 is always 0 - so we start the loop at index 1
      for k in range(1, numCoeffs):
       coeff = coeffs[k]
       if (coeff != 1):
         print(coeff, end='')
       if even_power:
         print('(2n + 1)', end='')
       self.__print_falling_factorial(k, 2*k)
       if even_power:
         print(f'/{(2*(2*k + 1))}', end='')
       else:
         print(f'/{2*k}', end='')

       if (k != numCoeffs - 1):
         print(' + ', end='')
    else:
      # Special case for 0
      print('(n + 1)', end='')
    print()

  def compute_sum_with_time_stat(self, power: int, n: int,
                                 stat: list[int]) -> mpz:
    sum: mpz = 0
    stat.clear()

    if (power < 0 or n < 0):
      stat.append(0)
      stat.append(0)
      return sum

    if power > 0:
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      coeffs: list[mpz] = self.__get_coefficients(power, n)
      after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      stat.append(after - before)
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      falling_factorial: mpz = 1
      even_power: bool = ((power & 1) == 0)
      numCoeffs = len(coeffs)
      # Coefficient 0 is always 0 - so we start the loop at index 1
      for k in range(1, numCoeffs):
        falling_factorial = falling_factorial*(n + k)*(n -k + 1)
        if (even_power):
          sum = sum + (coeffs[k]*falling_factorial*(2*n + 1))//(2*(2*k + 1))
        else:
          sum = sum + (coeffs[k]*falling_factorial)//(2*k)
    else:
      stat.append(0)
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      sum = n + 1

    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)
    return sum

   # The Central Factorial Numbers of second kind are as defined below:
   # T(2m, 2k) = k*k*T(2m - 2, 2k) + T(2m - 2, 2k - 2)
   # T(2m, 2m) = 1
   # We need only even coefficients.
  def __get_coefficients(self, power: int, max_n: int) -> list[mpz]:
    coeffs: list[mpz] = []

    if power < 0:
      return coeffs

    # Get half of power.  For odd power, round up half power
    m = (power >> 1) + (power & 1)
    i: int
    k: int
    T_2_2: mpz
    temp: mpz

    # We do not compute coefficients more than what is needed.  Example: power
    # is 1000 and the number of terms in the series is only 3 like
    # (1^1000 + 2^1000 + 3^1000).  In this case, falling factorials in many
    # terms will be 0.  There is no point in computing the coefficients for
    # these terms.

    max_num_coefficients = m + 1 # Array is 0 based, so add 1
    if (max_n < max_num_coefficients):
      max_num_coefficients = max_n + 1

    # Initialize the coefficient list
    for k in range(0, max_num_coefficients):
      coeffs.append(0)

    for i in range(0, m + 1):
      for k in range(0, max_num_coefficients):
        if (i == k):
          coeffs[k] = 1
        elif (i > 0 and k > 0):
          # Calculate T(2*i, 2*k) = T(2*i-2,2*k-2) + k*k*T(2*i-2, 2*k)
          # Use the previous values in coeffs[k] for all the computation before
          # setting the new value.
          temp = coeffs[k]*k*k + T_2_2
          T_2_2 = coeffs[k]
          coeffs[k] = temp
        else:
          coeffs[k] = 0
      T_2_2 = coeffs[0]
    return coeffs

  def __print_falling_factorial(self, start: int, num_terms: int) -> None:
    for i in range(num_terms):
      if (i == 0):
        print(f'(n + {start})', end='')
      else:
        diff = start - i
        if (diff == 0):
          print('n', end='')
        elif (diff > 0):
          print(f'(n + {diff})', end='')
        else:
          print(f'(n - {-diff})', end='')
