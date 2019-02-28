#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.qcdSmearProducer import *

mods = [
    qcdSmearProducer(),
]

files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_80X_v2_NanAOD_MC/PostProcess_v1/QCD_HT100to200/QCD_HT100to200_0.root"]
p=PostProcessor(".",files,cut=None, branchsel=None, haddFileName=True,outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_smear.txt",typeofprocess="smear",modules=mods,provenance=False)
p.run()
