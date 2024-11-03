theApp.EvtMax = -1

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = ['/afs/cern.ch/work/s/ssindhu/private/POWHEG_2024/events100k/test.pool.root']



from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD']

rivet.Analyses += [ 'tttt_event' ]
rivet.RunName = ''
rivet.HistoFile = 'tttt_event.yoda.gz' #output name
rivet.CrossSection = 11.91
#rivet.NominalWeightName = "_dyn__10_muR010000E+01_muF010000E+01_"
#rivet.IgnoreBeamCheck = True
#rivet.SkipWeights=True
job += rivet
