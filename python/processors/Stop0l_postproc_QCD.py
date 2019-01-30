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
from PhysicsTools.NanoSUSYTools.modules.tauMVAProducer import *

mods = [
    qcdSmearProducer(),
    #tauMVAProducer(),
]

files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]
#files=["root://cmseos.fnal.gov//store/user/benwu/Stop18/NtupleSyncMiniAOD/NanoSUSY/2018Xmas/prod2017MC_NANO.root"]
p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
p.run()
