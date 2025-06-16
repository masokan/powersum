#!/bin/sh
dir=`dirname $0`
scriptDir=`cd $dir; pwd`
cd $scriptDir
rm -f PowerSum.jar
javac com/github/masokan/powersum/*.java
jar -cfe PowerSum.jar com/github/masokan/powersum/PowerSumMain com/github/masokan/powersum/*.class
