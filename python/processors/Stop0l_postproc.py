#!/usr/bin/env python
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from PhysicsTools.NanoSUSYTools.modules.eleMiniCutIDProducer import *
#from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import *
#from PhysicsTools.NanoSUSYTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoSUSYTools.modules.DeepTopProducer import *
#from PhysicsTools.NanoSUSYTools.modules.updateGenWeight import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoSUSYTools.processors.ana_TrigNanoAOD import TrigEffAnalysis
#from PhysicsTools.NanoSUSYTools.processors.ana_TrigNanoAOD_2016 import TrigEffAnalysis_2016

def main(args):

    mods = [
        eleMiniCutID(),
        #Stop0lObjectsProducer(args.era),
        DeepTopProducer(args.era),
        #Stop0lBaselineProducer(args.era, args.isFastSim),
        # UpdateGenWeight(args.crossSection, args.nEvents)
	TrigEffAnalysis(args.era),
    ]

    #if args.era == "2016":
    #    mods += [
    #        puWeightProducer(pufile_mc,pufile_data,"pu_mc","pileup",verbose=False)
    #    ]

    #if args.era == "2017":
    #    mods += [
    #        puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)
    #    ]

    #files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]
    #files=["/eos/uscms/store/user/benwu/Stop18/NtupleSyncMiniAOD/NanoSUSY/2018Xmas/prod2017MC_NANO.root"]

    files = []
    lines = open(args.inputfile).readlines()
    for line in lines:
        files.append(line.strip())


    p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, modules=mods, histFileName="histos_METTrigNanoAOD.root",histDirName="metTrigAnalyzerMiniAOD")
    p.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-i', '--inputfile',
        default = "testing.txt",
        help = 'Path to the input filelist.')
    parser.add_argument('-o', '--outputfile',
                        default="./",
                        help = 'Path to the output file location.')
    parser.add_argument('-e', '--era',
        default = "", help = 'Year of production')
    parser.add_argument('-f', '--isFastSim', default = False)
    parser.add_argument('-d', '--isData', default = False)
    parser.add_argument('-c', '--crossSection')
    parser.add_argument('-n', '--nEvents', help = 'Number of Events')
    args = parser.parse_args()
    main(args)

