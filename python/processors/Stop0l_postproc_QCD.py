#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.eleMiniCutIDProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoSUSYTools.modules.qcdSmearProducer import *
#from PhysicsTools.NanoSUSYTools.modules.genjettest import *

#from PhysicsTools.NanoSUSYTools.modules.qcdskimmingfile import *
mods = [
    #qcdskimmingfile(),
      qcdSmearProducer(),
     # genjettest()
]

files=["../ttbar_v3_tree.root"]
#files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/HEMNano/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/2018_HEM_MC_RunIISpring18MiniAOD-100X_v10-v1-ext1/190105_203555/0000/BASE_80X_prodIso_NANO_143.root"]
#files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]

mods = [
    qcdSmearProducer(),
]

files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]
p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
p.run()
