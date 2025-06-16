#/bin/sh
# Script to run performance tests.  An optional argument specifies the language
# implementation to use to run the tests. The script runs tests with various
# powers and number of terms.  At the end, it plots the results and outputs
# encapsulated postscript files under the directory "perfTestResults-{language}
# The plots are for the total time (coefficient initialization + summation),
# summation time, and coefficient initialization time.
postprocess() {
  file="$1"
  prefix="$2"
  grep -E 'Computing S|Time taken' $file |
   sed -e 's/^Computing.*S.//' -e 's/[)]//' -e 's/Time taken = //' |
    awk -vprefix="$prefix" '
  BEGIN {
    ymins_total = "";
    ymaxs_total = "";
    ymin_total = 1145346789;
    ymax_total = 0;
    ymins_sum = "";
    ymaxs_sum = "";
    ymin_sum = 1145346789;
    ymax_sum = 0;
    faulhaberTotalFile = prefix "_faulhaber_total.csv";
    bernoulliTotalFile = prefix "_bernoulli_total.csv";
    stirlingTotalFile = prefix "_stirling_total.csv";
    eulerTotalFile = prefix "_euler_total.csv";
    centralTotalFile = prefix "_central_total.csv";
    faulhaberCoeffFile = prefix "_faulhaber_coeff.csv";
    bernoulliCoeffFile = prefix "_bernoulli_coeff.csv";
    stirlingCoeffFile = prefix "_stirling_coeff.csv";
    eulerCoeffFile = prefix "_euler_coeff.csv";
    centralCoeffFile = prefix "_central_coeff.csv";
    faulhaberSumFile = prefix "_faulhaber_sum.csv";
    bernoulliSumFile = prefix "_bernoulli_sum.csv";
    stirlingSumFile = prefix "_stirling_sum.csv";
    eulerSumFile = prefix "_euler_sum.csv";
    centralSumFile = prefix "_central_sum.csv";
    seriesFile = prefix "_series.csv";
    coeffTime["faulhaber"] = 0;
    coeffTime["bernoulli"] = 0;
    coeffTime["stirling"] = 0;
    coeffTime["euler"] = 0;
    coeffTime["central"] = 0;
  }
  {
    if ($0 ~ /.*,.*/) {
      numRecs = $0;
      sub("^.*,", "", numRecs);
      power = $0;
      sub(",.*$", "", power);
      printf("%s ", numRecs) >faulhaberTotalFile;
      printf("%s ", numRecs) >bernoulliTotalFile;
      printf("%s ", numRecs) >stirlingTotalFile;
      printf("%s ", numRecs) >eulerTotalFile;
      printf("%s ", numRecs) >centralTotalFile;
      printf("%s ", numRecs) >faulhaberSumFile;
      printf("%s ", numRecs) >bernoulliSumFile;
      printf("%s ", numRecs) >stirlingSumFile;
      printf("%s ", numRecs) >eulerSumFile;
      printf("%s ", numRecs) >centralSumFile;
      printf("%s ", numRecs) >seriesFile;
      count = 0;
    } else {
      numFields = split($0, t, ":");
      if (numFields < 3) {
        t[2] = t[1];
        t[3] = t[1];
      }
      switch(count) {
        case 0:
          print t[1] >faulhaberTotalFile;
          if (t[2] > coeffTime["faulhaber"]) {
            coeffTime["faulhaber"] = t[2];
          }
          print t[3] >faulhaberSumFile;
          break;
        case 1:
          print t[1] >bernoulliTotalFile;
          if (t[2] > coeffTime["bernoulli"]) {
            coeffTime["bernoulli"] = t[2];
          }
          print t[3] >bernoulliSumFile;
          break;
        case 2:
          print t[1] >stirlingTotalFile;
          if (t[2] > coeffTime["stirling"]) {
            coeffTime["stirling"] = t[2];
          }
          print t[3] >stirlingSumFile;
          break;
        case 3:
          print t[1] >eulerTotalFile;
          if (t[2] > coeffTime["euler"]) {
            coeffTime["euler"] = t[2];
          }
          print t[3] >eulerSumFile;
          break;
        case 4:
          print t[1] >centralTotalFile;
          if (t[2] > coeffTime["central"]) {
            coeffTime["central"] = t[2];
          }
          print t[3] >centralSumFile;
          break;
        case 5:
          print $0 >seriesFile;
          break;
      }
      val = t[1] + 0;
      if (val < ymin_total) {
        ymin_total = val;
        ymins_total = t[1];
      }
      if (val > ymax_total) {
        ymax_total = val;
        ymaxs_total = t[1];
      }
      val = t[3] + 0;
      if (val < ymin_sum) {
        ymin_sum = val;
        ymins_sum = t[3];
      }
      if (val > ymax_sum) {
        ymax_sum = val;
        ymaxs_sum = t[3];
      }
      count++;
    }

  }
  END {
    if (NR > 0) {
      printf("%d %d\n", power, int(coeffTime["faulhaber"])) >faulhaberCoeffFile;
      printf("%d %d\n", power, int(coeffTime["bernoulli"])) >bernoulliCoeffFile;
      printf("%d %d\n", power, int(coeffTime["stirling"])) >stirlingCoeffFile;
      printf("%d %d\n", power, int(coeffTime["euler"])) >eulerCoeffFile;
      printf("%d %d\n", power, int(coeffTime["central"])) >centralCoeffFile;
    }
    print ymins_total ":" ymaxs_total ";" ymins_sum ":" ymaxs_sum;
  }'
}

plotTimeVsN() {
  power=$1
  imageDir=$2
  plotDataDir=$3
  xrange=$4
  yrangeTotal=$5
  yrangeSum=$6
  plotFile=$plotDataDir/${power}_plotFileTotal
  cat >$plotFile <<EOT
set xlabel "Number of terms (n)"
set ylabel "Total CPU time (nanoseconds)"
set terminal png
set output '$imageDir/${power}_total.png'
set multiplot
set style data linespoints
set xrange [$xrange]
set yrange [$yrangeTotal]
set logscale x 10
set logscale y 10
set title "Power (m) = $power"
set arrow nohead lc rgb "#F748A5" from 10000000,960000 to 20000000,960000
set label "Faulhaber" at 22000000,960000
set arrow nohead lc rgb "#000000" from 10000000,480000 to 20000000,480000
set label "Bernoulli" at 22000000,480000
set arrow nohead lc rgb "#2271B2" from 10000000,240000 to 20000000,240000
set label "Stirling" at 22000000,240000
set arrow nohead lc rgb "#e69f00" from 10000000,120000 to 20000000,120000
set label "Euler" at 22000000,120000
set arrow nohead lc rgb "#359B73" from 10000000,60000 to 20000000,60000
set label "Central" at 22000000,60000
set arrow nohead lc rgb "#d55e00" from 10000000,30000 to 20000000,30000
set label "Series sum" at 22000000,30000
unset key
plot '$plotDataDir/${power}_faulhaber_total.csv' lc "#F748A5"
plot '$plotDataDir/${power}_bernoulli_total.csv' lc "#000000"
plot '$plotDataDir/${power}_stirling_total.csv' lc "#2271B2"
plot '$plotDataDir/${power}_euler_total.csv' lc "#e69f00"
plot '$plotDataDir/${power}_central_total.csv' lc "#359B73"
plot '$plotDataDir/${power}_series.csv' lc "#d55e00"
EOT
  gnuplot $plotFile

  plotFile=$plotDataDir/${power}_plotFileSum
  cat >$plotFile <<EOT
set xlabel "Number of terms (n)"
set ylabel "Sum CPU time (nanoseconds)"
set terminal png
set output '$imageDir/${power}_sum.png'
set multiplot
set style data linespoints
set xrange [$xrange]
set yrange [$yrangeSum]
set logscale x 10
set logscale y 10
set title "Power (m) = $power"
set arrow nohead lc rgb "#F748A5" from 10000000,960000 to 20000000,960000
set label "Faulhaber" at 22000000,960000
set arrow nohead lc rgb "#000000" from 10000000,480000 to 20000000,480000
set label "Bernoulli" at 22000000,480000
set arrow nohead lc rgb "#2271B2" from 10000000,240000 to 20000000,240000
set label "Stirling" at 22000000,240000
set arrow nohead lc rgb "#e69f00" from 10000000,120000 to 20000000,120000
set label "Euler" at 22000000,120000
set arrow nohead lc rgb "#359B73" from 10000000,60000 to 20000000,60000
set label "Central" at 22000000,60000
set arrow nohead lc rgb "#d55e00" from 10000000,30000 to 20000000,30000
set label "Series sum" at 22000000,30000
unset key
plot '$plotDataDir/${power}_faulhaber_sum.csv' lc "#F748A5"
plot '$plotDataDir/${power}_bernoulli_sum.csv' lc "#000000"
plot '$plotDataDir/${power}_stirling_sum.csv' lc "#2271B2"
plot '$plotDataDir/${power}_euler_sum.csv' lc "#e69f00"
plot '$plotDataDir/${power}_central_sum.csv' lc "#359B73"
plot '$plotDataDir/${power}_series.csv' lc "#d55e00"
EOT
  gnuplot $plotFile
}

printYMinMax() {
  awk '
    BEGIN {
      ymins = "";
      ymin = 32423134251;
      ymaxs = "";
      ymax = 0;
    }
    {
      if (FNR > 1) {
        y = $2 + 0;
        if (y < ymin) {
          ymin = y;
          ymins = $2;
        }
        if (y > ymax) {
          ymax = y;
          ymaxs = $2;
        }
      }
    }
    END {
     print ymins ":" ymaxs;
  }' "$@"
}

plotTimeVsM() {
  plotDataDir=$1
  imageDir=$2
  xrange=$3
  yrange=$4
  plotFile=$plotDataDir/plotFile_Power
  cat >$plotFile <<EOT
set xlabel "Power (m)"
set ylabel "Coeff. init. CPU time (nanoseconds)"
set terminal png
set output '$imageDir/coeffInitTime.png'
set multiplot
set style data linespoints
set xrange [$xrange]
set yrange [$yrange]
set logscale x 10
set logscale y 10
set title "Coeff. init. time versus power"
set arrow nohead lc rgb "#F748A5" from 5000,810000 to 10000,810000
set label "Faulhaber" at 12000,810000
set arrow nohead lc rgb "#000000" from 5000,270000 to 10000,270000
set label "Bernoulli" at 12000,270000
set arrow nohead lc rgb "#2271B2" from 5000,100000 to 10000,90000
set label "Stirling" at 12000,90000
set arrow nohead lc rgb "#e69f00" from 5000,30000 to 10000,30000
set label "Euler" at 12000,30000
set arrow nohead lc rgb "#359B73" from 5000,10000 to 10000,10000
set label "Central" at 12000,10000
unset key
plot '$plotDataDir/faulhaberCoeff.csv' lc "#F748A5"
plot '$plotDataDir/bernoulliCoeff.csv' lc "#000000"
plot '$plotDataDir/stirlingCoeff.csv' lc "#2271B2"
plot '$plotDataDir/eulerCoeff.csv' lc "#e69f00"
plot '$plotDataDir/centralCoeff.csv' lc "#359B73"
EOT
  gnuplot $plotFile
}

case $# in
  1) targetDir=$1
     lang=c
     ;;
  2) targetDir=$1
     lang=$2
     ;;
  *) echo "Usage: $0 <targetDir> [<language>]" 1>&2
     echo "<language> can be c or cpp; default is c"
     exit 1
     ;;
esac

scriptDir="`dirname $0`"
scriptDir="`cd $scriptDir; pwd`"

targetDir="`cd $targetDir; pwd`"

resultsDir=$targetDir/perfTestResults-${lang}
rm -rf $resultsDir
imageDir=$resultsDir/images
plotDataDir=$resultsDir/plotData
mkdir -p $imageDir $plotDataDir
srcDir=$scriptDir/../src/main/$lang
case $lang in
  c|cpp)  (cd $srcDir; make clean; make)
          if [ x$lang = "xc" ]
          then
              command=$srcDir/powersum
          else
              command=$srcDir/PowerSum
          fi
          ;;
  java)   (cd $srcDir; ./build.sh)
          command="java -jar $srcDir/PowerSum.jar"
          ;;
  python) rm -rf perftestenv
          python3 -m venv perftestenv
          . ./perftestenv/bin/activate
          pip3 install -r $srcDir/requirements.txt
          command="python3 $srcDir/power_sum_main.py"
          ;;
  *)      echo "Invalid language: $lang.  Should be one of c, cpp, java, or python" 1>&2
          exit 1
          ;;
esac

nrange="10:1000000000"
for power in 0010 0015 0020 0050 0100 0200 0500 1000 2000 3000 4000
do
  for terms in 10 25 50 75 100 250 500 750 1000 2000 5000 7500 10000 20000 \
               50000 75000 100000 200000 500000 750000 1000000
  do
     $command -sv $power $terms
  done  >${power}.log
  yrange=`postprocess ${power}.log $plotDataDir/$power`
  yrangeTotal=`echo $yrange|sed -e 's/;.*$//'`
  yrangeSum=`echo $yrange|sed -e 's/^.*;//'`
  plotTimeVsN $power $imageDir $plotDataDir $nrange $yrangeTotal $yrangeSum
  mv ${power}.log $plotDataDir
done

# Plot coeff init time versus power
mrange="10:50000"
cat $plotDataDir/[0-9]*_faulhaber_coeff.csv >>$plotDataDir/faulhaberCoeff.csv
cat $plotDataDir/[0-9]*_bernoulli_coeff.csv >>$plotDataDir/bernoulliCoeff.csv
cat $plotDataDir/[0-9]*_stirling_coeff.csv >>$plotDataDir/stirlingCoeff.csv
cat $plotDataDir/[0-9]*_euler_coeff.csv >>$plotDataDir/eulerCoeff.csv
cat $plotDataDir/[0-9]*_central_coeff.csv >>$plotDataDir/centralCoeff.csv
yrange=`printYMinMax $plotDataDir/*Coeff.csv`
plotTimeVsM $plotDataDir $imageDir $mrange $yrange

# Cleanup
echo "Cleaning up ..."
case $lang in
  c|cpp)  (cd $srcDir; make clean)
          ;;
  java)   (cd $srcDir; rm PowerSum.jar;
           find . -name '*.class'|xargs rm 2>/dev/null)
          command="java -jar $srcDir/PowerSum.jar"
          ;;
  python) rm -rf perftestenv
          (cd $srcDir; rm -rf __pycache__)
          ;;
esac
