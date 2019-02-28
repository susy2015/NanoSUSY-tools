#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.JetResSkim import *


mods = [
    JetResSkim(),
]

files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_80X_v2_NanAOD_MC/PostProcess_v1/QCD_HT100to200/QCD_HT100to200_0.root"]
p=PostProcessor(".",files,cut=None, branchsel=None,outputbranchsel="keep_and_drop_res.txt",typeofprocess="resp",modules=mods,provenance=False)
p.run()
