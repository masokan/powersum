#!/bin/sh
# This script runs functional tests on different language implementations. For
# each test, identical outputs are expected from all language implementations.
# For failed tests, the test logs are not deleted.  They can be used for
# debugging later.  The names of the log files have the unique test id embedded
# in them.

runOneTest() {
  testId=$1
  testOpts=$2
  outSuffix="_$testId.log"
  $cppSrcDir/PowerSum $testOpts 2>&1 |
    grep -E -v 'Time|Usage:|PowerSum' >cpp${outSuffix}
  $cSrcDir/powersum $testOpts 2>&1 |
    grep -E -v 'Time|Usage:|powersum' >c${outSuffix}
  java -jar $javaSrcDir/PowerSum.jar $testOpts 2>&1 |
    grep -E -v 'Time|Usage:|PowerSum.jar' >java${outSuffix}
  python3 $pythonSrcDir/power_sum_main.py $testOpts 2>&1 |
    grep -E -v 'Time|Usage:|power_sum_main.py' >python${outSuffix}

  cmp cpp${outSuffix} c${outSuffix} && rm c${outSuffix} &&
  cmp cpp${outSuffix} java${outSuffix} && rm java${outSuffix} &&
  cmp cpp${outSuffix} python${outSuffix} &&
    rm cpp${outSuffix} python${outSuffix}

  code=$?
  if [ $code -eq 0 ]
  then
    echo "Test Successful :-)"
  else
    echo "Test failed :-("
  fi
  return $code
}

scriptDir="`dirname $0`"
scriptDir="`cd $scriptDir; pwd`"
srcDir="$scriptDir/../src/main/"

cSrcDir=$srcDir/c
cppSrcDir=$srcDir/cpp
javaSrcDir=$srcDir/java
pythonSrcDir=$srcDir/python

# Print OS version
uname -srv

# Print processor information
cat /proc/cpuinfo

# Print compiler versions
gcc --version
g++ --version
java -version 2>&1
python3 -V

# Build the executables
(cd $cSrcDir; make clean; make)
(cd $cppSrcDir; make clean; make)
(cd $javaSrcDir; ./build.sh)

cd $scriptDir

# Create python virtual environment and install required packages
rm -rf functestenv
python3 -m venv functestenv
. ./functestenv/bin/activate
pip3 install -r $pythonSrcDir/requirements.txt

################################ TESTS START HERE #############################

# Test to verify coefficients
echo "Running test to verify coefficients for even power"
runOneTest 1 "-c 10"

echo "Running test to verify coefficients for 0th power"
runOneTest 2 "-c 0"

echo "Running test to verify coefficients for odd power"
runOneTest 3 "-c 17"

# Test to verify sum formula
echo "Running test to verify sum formula for odd power"
runOneTest 4 "-f 13"

echo "Running test to verify sum formula for 0th power"
runOneTest 5 "-f 0"

echo "Running test to verify sum formula for even power"
runOneTest 6 "-f 20"

# Test to verify sum
echo "Running test to verify sum for odd power"
runOneTest 7 "-s 23 100000"

echo "Running test to verify sum for even power"
runOneTest 8 "-s 26 1000"

echo "Running test to verify sum for 0th power"
runOneTest 9 "-s 0 1000"

echo "Running test to verify sum when number of terms is less than power"
runOneTest 10 "-s 19 7"

# Test to verify series sum
echo "Running test to verify sum with series sum up to large number"
runOneTest 11 "-sv 36 10000000"

echo "Running test to verify sum with series sum up to 1"
runOneTest 12 "-sv 13 1"

echo "Running test to verify sum with series sum up to 0 and non-zero power"
runOneTest 13 "-sv 22 0"

echo "Running test to verify sum with series sum up to 0 and power 0"
runOneTest 14 "-sv 0 0"

# Miscellaneous tests
echo "Running test with incorrect power"
runOneTest 15 "-sv -2 10"

echo "Running test with incorrect terms"
runOneTest 16 "-sv 15 -3877"

deactivate

# Cleanup
echo "Cleaning up ..."
rm -rf functestenv
(cd $cSrcDir; make clean)
(cd $cppSrcDir; make clean)
(cd $javaSrcDir; rm PowerSum.jar; find . -name '*.class'|xargs rm 2>/dev/null)
(cd $pythonSrcDir; rm -rf __pycache__)
