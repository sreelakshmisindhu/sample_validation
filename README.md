# Sample Validation



## Getting started
To clone the repository:
```
git clone ssh://git@gitlab.cern.ch:7999/ssindhu/sample_validation.git
```

## Generating grid pack

The folder `makegridpack` contains all the files for this. The jobOptions file is inside the `fourtops` folder. It can also be found here: [https://gitlab.cern.ch/atlas-physics/pmg/mcjoboptions/-/blob/master/950xxx/950512/mc.PhPy8EG_A14NNPDF23_fourtops_NLOQCDLOEW_valid.py](https://gitlab.cern.ch/atlas-physics/pmg/mcjoboptions/-/blob/master/950xxx/950512/mc.PhPy8EG_A14NNPDF23_fourtops_NLOQCDLOEW_valid.py)

To generate events without EW use the option `PowhegConfig.ewborn = 0`
Steps for setting up:
```
cd makegridpack
setupATLAS
asetup AthGeneration,23.6.39
Gen_tf.py --ecmEnergy=13000 --randomSeed=12345 --firstEvent=1 --maxEvents=10 --jobConfig=fourtops --outputTXTFile=test.lhe.tar.gz --outputEVNTFile=test.pool.root
```

This could take over 24 hours to run so it would be better to submit the job to HTCondor. The run and submit flies for this in also in the makegridpack folder. In the `run.sh` file change the path to the local path where you have cloned the repository to. submit the job using the command 

```
condor_submit submit.sub
```
Use ```condor_q``` to check the status of the jobs.

## Generating events
Events can be generated using the same `Gen_tf.py` command. To make sure the code uses the grid pack rename it as `mc_13TeV.PhPy8EG_A14NNPDF23_fourtops_NLOQCDLOEW_ath23.GRID.tar.gz` The part PhPy8EG_A14NNPDF23_fourtops_NLOQCDLOEW_ath23 can be changed according the sample but the format `mc_13TeV.[physicsShort].GRID.tar.gz` should be maintained.
Move this file to the fourtops folder. The grid pack has to be in the same folder as the JO for it to be identified. Run the Gen_tf command for as many events as necessary. This is the part that is being done centrally using the JIRA ticket: [https://its.cern.ch/jira/browse/ATLMCPROD-10204](https://its.cern.ch/jira/browse/ATLMCPROD-10204)

## Validation studies using Rivet

The Rivet routines are located in the rivet_analysis folder. This part is based on the tutorial for Rivet [https://topreco-tutorials.web.cern.ch/MCTutorial2023/RivetAndGeneration/ModifyJOsandRivet/](https://topreco-tutorials.web.cern.ch/MCTutorial2023/RivetAndGeneration/ModifyJOsandRivet/) Setting up Rivet:
```
cd makegridpack
setupATLAS
asetup AthGeneration,23.6.39
source setupRivet
```
To compile the routine run:
```
rivet-build tttt_event.so tttt_event.cc
```
To run Rivet using the Rivet routine use the `RunRivet_local_powheg.py` To test this locally with the samples produced on the grid download an example using the following command:
```
lsetup rucio
rucio download --nrandom 1 name_of_folder
```
Change the path for the input file in `RunRivet_local_powheg.py` to the path of the downloaded file. Run the plotting locally using:

```
athena RunRivet.py
```
To include the complete sample and run Rivet on the grid 
```
lsetup panda
pathena --extOutFile=MyOutput.yoda.gz --inDS=name_of_folder --outDS=user.`whoami`.origaS.v1 RunRivet_local_powheg.py
```
To compare with the MadGraph sample run Rivet on the Madgraph sample: `mc16_13TeV.412043.aMcAtNloPythia8EvtGen_A14NNPDF31_SM4topsNLO.merge.EVNT.e7101_e5984`
To test this locally download a random file from the file simalar to the Powheg sample and run Rivet on it.

To create plots use the command:
```

rivet-mkhtml --errs --no-weights  -o my_plots MyOutput.yoda.gz:"Title=powheg" MyOutput.yoda.gz:"Title=madgraph"
```


