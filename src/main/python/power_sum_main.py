from gmpy2 import mpz, mpq
from resource import *
import sys
import time
from central_factorial_power_sum import CentralFactorialPowerSum
from power_sum import PowerSum
from stirling_power_sum import StirlingPowerSum
from euler_power_sum import EulerPowerSum
from bernoulli_power_sum import BernoulliPowerSum
from faulhaber_power_sum import FaulhaberPowerSum

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

def print_usage(cmd: str):
  msg = '''Usage: {cmd} (-c|-f|-h|-s|-sv) [<power>] [<numTerms>]

<power> and <numTerms> should be greater than or equal to 0

Examples:
To print the help on usage:
{cmd}   -h or just   {cmd}

To print coefficients in the formula for power 10:
{cmd}   -c 10

To print formula for sum for power 5:
{cmd}   -f 5

To print the sum of series for power 6 for the first 20 terms:
The sum will be computed using the formula
{cmd}   -s 6 20

To compute the sum in two ways one using the formula and the
other with actual series expansion and verify the results for correctness
{cmd}   -sv 6 20\n
If <numTerms> is missing, a default of 20 is assumed'''.replace('{cmd}', cmd)
  eprint(msg)

def print_faulhaber_title():
  print('Faulhaber: -------------------------------------')

def print_bernoulli_title():
  print('Bernoulli: -------------------------------------')

def print_stirling_title():
  print('Stirling: --------------------------------------')

def print_euler_title():
  print('Euler: -----------------------------------------')

def print_central_factorial_title():
  print('Central Factorial: -----------------------------')

def print_cpu_time(before: int, after: int):
  print('Time taken = ', end='')
  cpuNano = after - before;
  print(cpuNano)

def print_coefficients(coeffs: list[mpq]) -> None:
  for i in range(len(coeffs)):
    print(f' {coeffs[i]}', end='')
  print('')

def get_coefficients_timed(ps: PowerSum, power: int):
  before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
  coeffs = ps.get_coefficients(power)
  after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
  print_coefficients(coeffs)
  print_cpu_time(before, after)

def compute_and_print_sum_timed(ps: PowerSum, power: int,
                                num_terms: int) -> mpz:
  stat: list[int] = []
  sum: mpz = ps.compute_sum_with_time_stat(power, num_terms, stat)
  print(f'Sum computed = {sum}')
  total_time = stat[0] + stat[1]
  print('Time taken = ', end='')
  print(f'{total_time}:{stat[0]}:{stat[1]}')
  return sum

args = sys.argv
cmd = 'python ' + args[0]
if len(args) < 2 or args[1] == '-h':
  print_usage(cmd)
  sys.exit(0)

valid_options = ['-c', '-f', '-s', '-sv']
power = 1
num_terms = 20

if len(args) > 2:
  power = int(args[2])
  if power < 0:
    eprint('Invalid power: ' + args[2])
    print_usage(cmd)
    sys.exit(1)

if len(args) > 3:
  num_terms = int(args[3])
  if num_terms < 0:
    eprint('Invalid number of terms: ' + args[3])
    print_usage(cmd)
    sys.exit(1)

if args[1] in valid_options:
  # Set the highest limit for printing sum values which can be thousands of
  # digits
  sys.set_int_max_str_digits(0)
  fps = FaulhaberPowerSum()
  bps = BernoulliPowerSum()
  sps = StirlingPowerSum()
  eps = EulerPowerSum()
  cps = CentralFactorialPowerSum()
  match args[1]:
    case '-c':
      print(f'Computing coefficients for power {power}')
      print_faulhaber_title()
      get_coefficients_timed(fps, power)
      print_bernoulli_title()
      get_coefficients_timed(bps, power)
      print_stirling_title()
      get_coefficients_timed(sps, power)
      print_euler_title()
      get_coefficients_timed(eps, power)
      print_central_factorial_title()
      get_coefficients_timed(cps, power)
    case '-f':
      print_faulhaber_title()
      fps.print_sum_formula(power)
      print_bernoulli_title()
      bps.print_sum_formula(power)
      print_stirling_title()
      sps.print_sum_formula(power)
      print_euler_title()
      eps.print_sum_formula(power)
      print_central_factorial_title()
      cps.print_sum_formula(power)
    case '-s'|'-sv':
      print(f'Computing S({power}, {num_terms})');
      print_faulhaber_title()
      faulhaber_sum: mpz = compute_and_print_sum_timed(fps, power, num_terms)
      print_bernoulli_title()
      bernoulli_sum: mpz = compute_and_print_sum_timed(bps, power, num_terms)
      print_stirling_title()
      stirling_sum: mpz = compute_and_print_sum_timed(sps, power, num_terms)
      print_euler_title()
      euler_sum: mpz = compute_and_print_sum_timed(eps, power, num_terms)
      print_central_factorial_title()
      cf_sum: mpz = compute_and_print_sum_timed(cps, power, num_terms)
      if (args[1] == '-sv'):
        print('Series addition:--------------------------------')
        before = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
        sum_from_series:mpz = sps.compute_sum_using_series(power, num_terms)
        after = time.clock_gettime_ns(time.CLOCK_THREAD_CPUTIME_ID)
        print(f'Sum computed = {sum_from_series}')
        print_cpu_time(before, after)
        if faulhaber_sum == sum_from_series:
          print('The sum matches with Faulhaber formula :-)')
        else:
          print('The sums do not match for Faulhaber formula :-(')
        if bernoulli_sum == sum_from_series:
          print('The sum matches with Bernoulli formula :-)')
        else:
          print('The sums do not match for Bernoulli formula :-(')
        if stirling_sum == sum_from_series:
          print('The sum matches with Stirling formula :-)')
        else:
          print('The sums do not match for Stirling formula :-(')
        if euler_sum == sum_from_series:
          print('The sum matches with Euler formula :-)')
        else:
          print('The sums do not match for Euler formula :-(')
        if cf_sum == sum_from_series:
          print('The sum matches with Central Factorial formula :-)')
        else:
          print('The sums do not match for Central Factorial formula :-(')
  sys.exit(0)
else:
  print_usage(cmd)
  sys.exit(1)
