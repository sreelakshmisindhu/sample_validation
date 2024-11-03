# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

#--------------------------------------------------------------
# This is an example joboption to generate events with Powheg
# using ATLAS' interface. Users should optimise and carefully
# validate the settings before making an official sample request.
#--------------------------------------------------------------

#--------------------------------------------------------------
# EVGEN configuration
#--------------------------------------------------------------
evgenConfig.description = "POWHEG+Pythia8 4 tops production with A14 NNPDF2.3 tune."
evgenConfig.keywords = ["SM", "top"]
evgenConfig.contact = ["tpelzer@cern.ch"]

# --------------------------------------------------------------
# Load ATLAS defaults for the Powheg 4 tops process
# --------------------------------------------------------------
include("PowhegControl/PowhegControl_fourtops_Common.py")
PowhegConfig.PDF = 303400 
PowhegConfig.mu_R = 1
PowhegConfig.mu_F = 1
PowhegConfig.btlscalereal = 1
PowhegConfig.btlscalect = 1
PowhegConfig.withdamp = 1
PowhegConfig.bornzerodamp = 1
PowhegConfig.dynamic_hdamp = 1

# --------------------------------------------------------------
# Settings
# --------------------------------------------------------------
# define the decay mode
PowhegConfig.decay_mode = "t t~ t t~ > all" # inclusive is the default
#PowhegConfig.decay_mode = "t t~ t t~ > 4l" # exactly 4 leptons
#PowhegConfig.decay_mode = "t t~ t t~ > 3l4l" # at least 3 leptons
#PowhegConfig.decay_mode = "t t~ t t~ > 2lSS" # exactly 2 leptons, same sign
#PowhegConfig.decay_mode = "t t~ t t~ > 2lOS" # exactly 2 leptons, opposite sign
#PowhegConfig.decay_mode = "t t~ t t~ > 1l" # exactly 1 lepton
#PowhegConfig.decay_mode = "t t~ t t~ > allhad" # fully-hadronic
#PowhegConfig.decay_mode = "t t~ t t~ > undecayed" # undecayed (e.g. to produce a lhe file)
## for handling decays with MadSpin
#PowhegConfig.decay_mode = "t t~ t t~ > all [MadSpin]"
#PowhegConfig.MadSpin_decays= ["decay t > w+ b, w+ > all all", "decay t~ > w- b~, w- > all all"]
#PowhegConfig.MadSpin_process= "generate p p > t t~ t t~ [QCD]" # this process is default - can be changed (for studies)
## additional MadSpin parameters to improve the integration
#PowhegConfig.MadSpin_max_weight_ps_point= 1000
#PowhegConfig.MadSpin_Nevents_for_max_weight= 250

# --------------------------------------------------------------
# Generate events
# --------------------------------------------------------------
PowhegConfig.ncall1       = "16000"
PowhegConfig.ncall2       = "12000"
PowhegConfig.ncall1rm     = "4000"
PowhegConfig.ncall2rm     = "4000"
PowhegConfig.nubound      = "8000"
PowhegConfig.generate()

#--------------------------------------------------------------
# Pythia8 showering with the A14 NNPDF2.3 tune, main31 routine
#--------------------------------------------------------------
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("Pythia8_i/Pythia8_Powheg_Main31.py")
# Setting the appropriate number of final state particles for the main31 routine
genSeq.Pythia8.Commands += [ 'Powheg:NFinal = 4' ]
