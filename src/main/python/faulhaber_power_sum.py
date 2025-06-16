from gmpy2 import mpz, mpq
from power_sum import PowerSum
import time

class FaulhaberPowerSum(PowerSum):

  def __init__(self) -> None:
    return

  # There is no simple recurrence relation to generate the coefficients.
  # According to A. W. F. Edwards
  # (http://www.pietrocola.eu/Fontecchio2019/A%20quick%20route%20to%20sums%20of%20powers%20by%20A.W.F.Edwards%20(1).pdf),
  # they can be obtained by matrix inversion where matrix rows are initialized
  # using the methods outlined in the paper.  In general, a matrix requires
  # O(m^2) storage and the inversion requires O(m^p) time where p > 2.
  # Instead of matrix inversion, we consider the problem as solving a system of
  # linear equtaions.  So, we contruct an augmented matrix and transform it to
  # reduced row echelon (rre) form.  In our case, we notice that if we properly
  # construct the matrix, it can be a row echelon matrix.
  # All we need to do is transform it to a rre matrix and solve for the
  # coefficients.  The transformation can be achieved sequentially one row at a
  # time.  There is no need to hold the entire matrix in memory thus reducing
  # the storage requirement to O(m).  The time complexity is also brought down
  # to O(m^2).
  def get_coefficients(self, power: int) -> list[mpq]:
    coeffs: list[mpq] = []
    if power < 0:
      return coeffs

    first_row: list[mpz] = []
    next_row: list[mpz] = []
    n_limit: int
    pivot_index: int
    odd_power: bool = False
    scale_by: mpz
    next_scale_by: mpz

    if ((power & 1) == 1):
      odd_power = True
      n_limit = (power + 1)//2
    else:
      n_limit = power//2 + 1

    # Create the first row.  Each row is created with an augmented column entry
    # of 1
    if (odd_power):
      scale_by = self.__create_row_for_odd_power(n_limit, n_limit, first_row)
    else:
      scale_by = self.__create_row_for_even_power(n_limit, n_limit, first_row)

    coeffs.append(mpq(first_row[-1], scale_by))

    pivot_index = 1
    # Create subsequent rows and generate one coefficient at a time
    for i in range(n_limit - 1, 0, -1):
      # Reinitialize the augmented column
      first_row[-1] = 0
      pivot: mpz = first_row[pivot_index]
      if (pivot == 0):
        # Once we reach a 0, all other columns will be zero.  This means all
        # other coeffs will be zero
        break
      next_row.clear()
      if (odd_power):
        next_scale_by = self.__create_row_for_odd_power(n_limit, i, next_row)
      else:
        next_scale_by = self.__create_row_for_even_power(n_limit, i, next_row)
      scale_by *= next_scale_by

      # Get the value in the pivot column in the next row
      pivot_in_next: mpz = next_row[pivot_index]

      # Multiply rest of the columns that follow the pivot in the first row by
      # next pivot column and the next row by pivot and subtract from first
      for j in range(pivot_index + 1, len(next_row)):
        first_row[j] *= pivot_in_next
        first_row[j] -= pivot*next_row[j]
      coeffs.append(mpq(first_row[-1], scale_by))
      pivot_index += 1

    return coeffs

  def print_sum_formula(self, power: int) -> None:
    if power > 0:
      coeffs: list[mpq] = self.get_coefficients(power)

      if ((power & 1) == 0):
        print('(2n + 1)', end='')
      exponent = (power + 1)//2
      print('{', end='')
      for i in range(0, len(coeffs)):
        if (i != 0):
          print(' + ', end='')
        print(f'({coeffs[i]})N', end='')
        if (exponent != 1):
          print(f'^{exponent}', end='')
        exponent -= 1
      print('}/2')
      print('where N = n(n + 1)')
    else:
      if power == 0:
        print('(n + 1)')

  def compute_sum_with_time_stat(self, power: int, n: int,
                                 stat: list[int]) -> mpz:
    result: mpz = 0
    stat.clear()
    if (power < 0 or n < 0):
      stat.append(0)
      stat.append(0)
      return result

    if power > 0:
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      coeffs: list[mpq] = self.get_coefficients(power)
      after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      stat.append(after - before)
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      if (n > 0):
        sum: mpq = 0
        N: mpz = n
        N *= (n + 1)
        NPow: mpz = N

        # For odd powers, the first term starts with N^2.  There is no N term
        if ((power & 1) == 1):
          NPow *= N
        for i in range(len(coeffs) - 1, -1, -1):
          sum += coeffs[i]*NPow
          NPow *= N
        if ((power & 1) == 0):
          sum = sum*(2*n + 1)
        (result, denominator) = sum.as_integer_ratio()
        result //= 2
    else:
      stat.append(0)
      before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
      result = n + 1
    after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
    stat.append(after - before)
    return result

  # Reverse the row and return the first non-zero column value
  def __post_process_row(self, row: list[mpz]) -> mpz:
    row.reverse()
    # Augment the row with 1 which corresponds to the coefficient of the
    # power of N (Note: N = n(n+1))
    row.append(1)
    scale_by: mpz = 1
    for j in range(0, len(row)):
      if (row[j] != 0):
        scale_by = row[j]
        break
    return scale_by

  def __create_row_for_even_power(self, n_limit: int, row_num: int,
                                  row: list[mpq]) -> mpz:
    num_added: int = 0
    ncr1: mpz = 0
    ncr2: mpz = 0
    for i in range(2*row_num - 1, 0, -2):
      if (i > row_num):
        row.append(0)
      else:
        # Calculate nCr only once and incrementally update it in the loop to
        # avoid redundant computation that will increase the time by an order
        # of magnitude
        if (ncr1 == 0):
          ncr1 = self.nCr(row_num, i)
        else:
          ncr1 *= (i + 2)*(i + 1)
          ncr1 //= (row_num - i - 1)*(row_num - i)
        if (ncr2 == 0):
          ncr2 = self.nCr(row_num - 1, i)
        else:
          ncr2 *= (i + 2)*(i + 1)
          ncr2 //= (row_num - i - 2)*(row_num - i - 1)
        row.append(ncr1 + ncr2)
      num_added += 1

    while(num_added < n_limit):
      row.append(0)
      num_added += 1

    return self.__post_process_row(row)

  def __create_row_for_odd_power(self, n_limit: int, row_num: int,
                                 row: list[mpq]) -> mpz:
    num_added: int = 0
    ncr: mpz = 0
    for i in range(2*row_num - 1, 0, -2):
      if (i > row_num):
        row.append(0)
      else:
        # Calculate nCr only once and incrementally update it in the loop to
        # avoid redundant computation that will increase the time by an order
        # of magnitude
        if (ncr == 0):
          ncr = self.nCr(row_num, i)
        else:
          ncr *= (i + 2)*(i + 1)
          ncr //= (row_num - i - 1)*(row_num - i)
        row.append(ncr)
      num_added += 1

    while(num_added < n_limit):
      row.append(0)
      num_added += 1

    return self.__post_process_row(row)
