#!/bin/bash

runAllTemplate=runAll_oneSim.sh
# nonIntroFracs=$(echo .2 .4 .6 .8)
# introTimes=$(echo .1 .5 .75)
headDir=$PWD

# for frac in $nonIntroFracs
# do
#   for time in $introTimes
#   do
#     introProb=$(bc <<< "scale=2; 1-$frac")
#     base=msSims_$time\_1000000_$introProb
#     mkdir $base
#     sed "s:NONINTROFRAC:$frac:g" $runAllTemplate | \
#     sed "s:INTROTIME:$time:g" | sed "s:BASENAME:$base:g" > $base/runAll.$base.sh
#     cd $base
#     mkdir code
#     ln -s $headDir/code/* code
#     chmod u+x runAll.$base.sh
#     chmod u+x code/*
#     sbatch ./runAll.$base.sh
#     cd $headDir
#   done
# done

overallOut=msSims.allSimulations.results
echo "time,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > msSims.allSimulations.results

###some didn't finish
for run in .5,.2 .75,.2 .75,.4
do
  time=$(echo $run|cut -d "," -f 1)
  introProb=$(echo $run|cut -d "," -f 2)
  base=msSims_$time\_1000000_$introProb
  sed "s:NONINTROFRAC:$frac:g" $runAllTemplate | \
  sed "s:INTROTIME:$time:g" | sed "s:BASENAME:$base:g" > $base/runAll.$base.sh
  cd $base
  mkdir code
  ln -s $headDir/code/* code
  chmod u+x runAll.$base.sh
  chmod u+x code/*
  sbatch ./runAll.$base.sh
  cd $headDir
done


#run after the rest is done
# for test in msSims_*_*
# do cd $test
# fdOut=$test.popStats.boots.summary
# quiblOut=$test.combined.quibl.out
# ../code/collateResults.slurm $test $fdOut $quiblOut $headDir/$overallOut
# cd $headDir
# done
