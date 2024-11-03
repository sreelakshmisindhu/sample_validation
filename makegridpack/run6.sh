#!/bin/bash

echo "welcome to HTCondor"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export AtlasSetup=/afs/cern.ch/atlas/software/dist/AtlasSetup
#alias asetup='source $AtlasSetup/scripts/asetup.sh'

#export ATHENA_PROC_NUMBER=10

source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
source $AtlasSetup/scripts/asetup.sh AthGeneration,23.6.38,here
#pkgco.py -A MadGraphControl
#cd Generators/MadGraphControl/cmt
#make clean; make
cd $TestArea
mkdir run6/
cd run6/

cp -r /afs/cern.ch/work/s/ssindhu/private/POWHEG_2024/October2024_tests/fourtops .
#cp /afs/cern.ch/user/p/peiffer/work/MadgraphTEST/run/$1 .

ls -l

echo "--------- START ---------"

Gen_tf.py --ecmEnergy=13000 --randomSeed=12345 --firstEvent=1 --maxEvents=10 --jobConfig=fourtops --outputTXTFile=test.lhe.tar.gz --outputEVNTFile=test.pool.root

echo "--------- END ---------"
cd ..
cp -r run6 /afs/cern.ch/work/s/ssindhu/private/POWHEG_2024/October2024_tests/
#cp *events.tar.gz /afs/cern.ch/work/s/ssindhu/private/POWHEG_2024/October2024_tests/

#cp log.generate /afs/cern.ch/user/p/peiffer/.
