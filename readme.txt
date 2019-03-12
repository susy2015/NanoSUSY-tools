Smearing QCD Notes NANOAOD
In PhysicsTools/NanoSUSYTools/python/processors/
You need to create a QCD file with all of the qcd_orig files that you want to run over.
An example of all of these files is located in my area: /uscms_data/d3/mkilpatr/CMSSW_10_2_9/src/PhysicsTools/NanoSUSYTools/python/processors

> python Stop0l_postproc_res.py

For Condor:
> python SubmitLPC.py -f ../Stop0l_postproc_res.py -c ../../../../../StopCfg/sampleSets_postProcess_2016_QCD.cfg -o /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod

hadd the output files together with the name "jetResSkim_combined_filtered_CHEF_NANO.root"

> mv jetResSkim_combined_filtered_CHEF_NANO.root /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod/.

You need to get the jet Res for the Tail Smear
> root -l -b -q ../rootlogon.C GetJetResForTailSmear.C+

> mkdir skims
> cd skims
> cp ../JetResDiagnostic.C .

Use JetResDiagnostic.C to create the pngs from the previous command ^
> root -l -b -q ../../rootlogon.C JetResDiagnostic.C+

> python Stop0l_postproc_QCD.py

For Condor:
> python SubmitLPC.py -f ../Stop0l_postproc_QCD.py -c ../../../../../StopCfg/sampleSets_preProcess_2016_QCD.cfg -o /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod/

After QCD smearing you need to run the add weight part of the NTuples.
You need to now create trees to run over and create the SF for the files. 
Once you create the trees from this file you need to have the met_tree.root, ttbarplusw_tree.root, and qcd_tree.root in the same directory.
You need to link the directory to the input files in the MakeQCDRespTailSF.C in the AnalysisMethods/macros/JetMETStudies/ directory.

> cd AnalysisMethods/macros/JetMETStudies/
> root -l -b -q ../rootlogon.C MakeQCDRespTailSF.C+

This will run over the input files and calculate the correct SF for the QCD. Once it finishes there will be an output root file. You need to put this file in the data directory under corrections/ and in the right year. 

After getting the SF the QCD needs to be postprocessed one last time.
