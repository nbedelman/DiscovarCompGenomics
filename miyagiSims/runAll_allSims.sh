#!/bin/bash

# master script to run analysis of all simulation sets

### setup global environment and define global variables
module load centos6/0.0.1-fasrc01  phyml/2014Oct16-fasrc01
module load Anaconda/5.0.1-fasrc01
source activate quibl


phymlPath=/n/sw/fasrcsw/apps/Core/phyml/2014Oct16-fasrc01/bin/phyml

headDir=$PWD
dataDir=$headDir/data

quiblTemplate=quiblInput.template.txt

### run code
##ancient_tests
#
# mkdir ancient_tests
# ancientData=ancient_tests.allSimulations.results
# windSize=1000
# echo "time,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $ancientData
#
# for time in $(ls $dataDir/revision_dem/ancient_tests/|grep -v seq-gen)
# do
#   for popSize in $(ls $dataDir/revision_dem/ancient_tests/$time/sgOut)
#   do
#   for frac in $(ls $dataDir/revision_dem/ancient_tests/$time/sgOut/$popSize)
#   do
#     testDirectory=ancient_tests/$time\_$popSize\_$frac
#     mkdir $testDirectory
#
#     cp $quiblTemplate $testDirectory
#     base=$time\_$popSize\_$frac
#     testName=$time
#     testDataDir=$dataDir/revision_dem/ancient_tests/$time/sgOut/$popSize/$frac
#
#     sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
#     sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$ancientData:" | \
#     sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:"> $testDirectory/runAll.$base.sh
#
#     cd $testDirectory
#     chmod u+x runAll.$base.sh
#     ./runAll.$base.sh
#     cd $headDir
# done
# done
# done

##ancient_repeat
#
# ancientData=ancient_tests.allSimulations.results
# windSize=1000
# echo "time,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $ancientData
# #
# for time in $(ls $dataDir/ancient_tests/sgRepeat|grep -v seq-gen)
# do
#   for popSize in $(ls $dataDir/ancient_tests/sgRepeat/$time)
#   do
#   for frac in $(ls $dataDir/ancient_tests/sgRepeat/$time/$popSize)
#   do
#     testDirectory=ancient_tests/$time\_$popSize\_$frac
#     mkdir $testDirectory
#
#     cp $quiblTemplate $testDirectory
#     base=$time\_$popSize\_$frac
#     testName=$time
#     testDataDir=$dataDir/ancient_tests/sgRepeat/$time/$popSize/$frac
#
#     sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
#     sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$ancientData:" | \
#     sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:"> $testDirectory/runAll.$base.sh
#
#     cd $testDirectory
#     chmod u+x runAll.$base.sh
#     ./runAll.$base.sh
#     cd $headDir
# done
# done
# done
#

##demo_tests
#
# mkdir demo_tests
# demoData=demo_tests.allSimulations.results
# windSize=1000
# echo "test,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $demoData
#
# #these are not all named uniformly, so just make a list
# demos=$(echo int_dec_anc  int_inc_anc  null rec_dec  rec_inc  src_dec	src_inc)
#
# for demo in $demos
# do
#   for popSize in $(ls $dataDir/revision_dem/demography_tests/$demo/sgOut)
#   do
#   for frac in $(ls $dataDir/revision_dem/demography_tests/$demo/sgOut/$popSize)
#   do
#     testDirectory=demo_tests/$demo\_$popSize\_$frac
#     mkdir $testDirectory
#
#     cp $quiblTemplate $testDirectory
#     base=$demo\_$popSize\_$frac
#     testName=$demo
#     testDataDir=$dataDir/revision_dem/demography_tests/$demo/sgOut/$popSize/$frac
#
#     sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
#     sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$demoData:" | \
#     sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:" > $testDirectory/runAll.$base.sh
#
#     cd $testDirectory
#     chmod u+x runAll.$base.sh
#     ./runAll.$base.sh
#     cd $headDir
# done
# done
# done

##demoData_redo (250kb Ne)
# demoData=demo_tests.allSimulations.results
# windSize=1000
# echo "test,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $demoData
#
# #these are not all named uniformly, so just make a list
# demos=$(echo int_dec_anc  int_inc_anc  null rec_dec  rec_inc  src_dec	src_inc)
#
# for demo in $demos
# do
#   for popSize in $(ls $dataDir/demography_tests/sgRepeat_demo/$demo)
#   do
#   for frac in $(ls $dataDir/demography_tests/sgRepeat_demo/$demo/$popSize)
#   do
#     testDirectory=demo_tests/$demo\_$popSize\_$frac
#     mkdir $testDirectory
#
#     cp $quiblTemplate $testDirectory
#     base=$demo\_$popSize\_$frac
#     testName=$demo
#     testDataDir=$dataDir/demography_tests/sgRepeat_demo/$demo/$popSize/$frac
#
#     sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
#     sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$demoData:" | \
#     sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:" > $testDirectory/runAll.$base.sh
#
#     cd $testDirectory
#     chmod u+x runAll.$base.sh
#     ./runAll.$base.sh
#     cd $headDir
# done
# done
# done
#

##5KB ancient_tests

# mkdir -p 5KB_ancient_tests
summaryFile=5KB_ancient_tests.allSimulations.results
windSize=5000
echo "test,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $summaryFile

for time in $(ls $dataDir/5KB_anc/|grep -v seq-gen)
do
  for popSize in $(ls $dataDir/5KB_anc/$time/sgOut)
  do
  for frac in $(ls $dataDir/5KB_anc/$time/sgOut/$popSize)
  do
    testDirectory=5KB_ancient_tests/$time\_$popSize\_$frac
    if [ ! -d $testDirectory ];then
    mkdir -p $testDirectory

    cp $quiblTemplate $testDirectory
    base=$time\_$popSize\_$frac
    testName=$time
    testDataDir=$dataDir/5KB_anc/$time/sgOut/$popSize/$frac

    sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
    sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$summaryFile:" | \
    sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:"> $testDirectory/runAll.$base.sh

    cd $testDirectory
    chmod u+x runAll.$base.sh
    ./runAll.$base.sh &
    cd $headDir
  fi
done
done
done


####recomb tests
# mkdir -p recomb_tests
# summaryFile=recomb_tests.allSimulations.results
# windSize=1000
# echo "test,popSize,frac,fd,fdMax,fdMin,D,DMax,DMin,C,maxC,minC,lambda,maxLambda,minLambda,introProp,maxIntroProp,minIntroProp" > $summaryFile
# # #
# # #for test in recomb_5 recomb_10
# for test in recomb_11
# do
#   for popSize in $(ls $dataDir/$test/sgOut)
#   do
#   for frac in $(ls $dataDir/$test/sgOut/$popSize)
#   do
#     testDirectory=recomb_tests/$test\_$popSize\_$frac
#     mkdir -p $testDirectory
#
#     cp $quiblTemplate $testDirectory
#     base=$test\_$popSize\_$frac
#     testName=$test
#     testDataDir=$dataDir/$test/sgOut/$popSize/$frac
#
#     sed "s:TESTDATADIR:$testDataDir:" runAll_oneTrial_template.sh  | \
#     sed "s:BASENAME:$base:" | sed "s:OVERALLOUT:$headDir/$summaryFile:" | \
#     sed "s:HEADDIR:$headDir:" |sed "s:WINDSIZE:$windSize:"> $testDirectory/runAll.$base.sh
#
#     cd $testDirectory
#     chmod u+x runAll.$base.sh
#     ./runAll.$base.sh &
#     cd $headDir
# done
# done
# done
