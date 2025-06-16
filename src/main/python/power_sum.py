from abc import ABC, abstractmethod
from gmpy2 import mpz, mpq, powmod

# PowerSum base class
class PowerSum(ABC):

  @abstractmethod
  # To get the coefficients in the sum formula.
  # Parameters:
  #   power - desired power (IN)
  # Return value
  #   a list of coefficients
  #
  def get_coefficients(self, power: int) -> list[mpq]:
    pass

  @abstractmethod
  # To print the sum formula in stdout for a given power
  # Parameters:
  #   power - desired power (IN)
  #
  def print_sum_formula(self, power: int) -> None:
    pass

  @abstractmethod
  # To compute the sum for a specific power and the number of terms and obtain
  # CPU times
  # Parameters:
  #   power - desired power (IN)
  #   n - number of terms (IN)
  #   stat - list in which the coefficient initialization time and
  #          summation time (in nanosec.) are returned (IN/OUT)
  # Return value
  #   sum computed
  #
  def compute_sum_with_time_stat(self, power: int, n: int,
                                 stat: list[int]) -> mpz:
    pass

  # To compute the sum for a specific power and the number of terms
  # Parameters:
  #   power - desired power (IN)
  #   n - number of terms (IN)
  # Return value
  #   sum computed
  #
  def compute_sum(self, power: int, n: int) -> mpz:
    stat: list[int] = []
    return self.compute_sum_with_time_stat(power, n, stat)

  def to_rational(self, coeffs: list[mpz]) -> list[mpq]:
    out_coeffs: list[mpq] = [mpq(c, 1) for c in coeffs]
    return out_coeffs

  def compute_sum_using_series(self, power: int, n: int) -> mpz:
    sum: mpz = 0
    if (power < 0 or n < 0):
      return sum
    for count in range(0, n + 1):
      result: mpz = pow(count, power)
      sum += result
    return sum

  def nCr(self, n: int, r: int) -> mpz:
    num: int = n

    if (r > n):
      return 0
    if (r > n/2):
      r = n - r
    if (r == 0):
      return 1

    result: mpz = num
    for i in range(2, r + 1):
      num -= 1
      result *= num
      result //= i
    return result
