#!/usr/bin/env python2.7
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gROOT,gSystem,gDirectory
from scharmfit.utils import load_susyfit
load_susyfit()
from ROOT import ConfigMgr,FitConfig #this module comes from gSystem.Load("libSusyFitter.so")
gROOT.Reset()

from ROOT import TFile, RooWorkspace, TObject, TString, RooAbsReal, RooRealVar, RooFitResult, RooDataSet, RooAddition, RooArgSet,RooAbsData,RooRandom 
from ROOT import Util, TMath
from ROOT import RooFit
from ROOT import RooExpandedFitResult
    
import os
import sys
from sys import exit

from SysTableTex import *
import pickle

from logger import Logger
log = Logger('SysTable')
#log.setLevel('DEBUG')

def getnamemap():

  namemap = {}

  return namemap

def _slim_regions(regionList, regionCat):
  """
  remove the regions that don't exits.
  God do I hate HistFitter... I can't vouch for the correctness
  of anything that comes out of it (sorry LHC).
  """
  keep_regions = []
  for region in regionList:
    full_name = Util.GetFullRegionName(regionCat, region)
    if str(full_name):
      keep_regions.append(region)
  return keep_regions
  

def latexfitresults( filename, region='3jL', sample='', resultName="RooExpandedFitResult_afterFit", dataname='obsData', doAsym=True):

  namemap = {} ## add this if I want description
  namemap = getnamemap() ## add this if I want description

  ############################################
  workspacename = 'w'
  w = Util.GetWorkspaceFromFile(filename,workspacename)

  if w==None:
    print "ERROR : Cannot open workspace : ", workspacename
    sys.exit(1) 

  result = w.obj(resultName)
  if result==None:
    print "ERROR : Cannot open fit result ", resultName
    sys.exit(1)

  snapshot =  'snapshot_paramsVals_' + resultName
  w.loadSnapshot(snapshot)

  data_set = w.data(dataname)
  if data_set==None:
    print "ERROR : Cannot open dataset : ", "data_set"
    sys.exit(1)
      
  regionCat = w.obj("channelCat")
  data_set.table(regionCat).Print("v");

  regionFullName = Util.GetFullRegionName(regionCat, region);

  chosenSample = False
  if sample is not '':
    chosenSample = True
        
  #####################################################

  regSys = {}

  regionCatStr = 'channelCat==channelCat::' + regionFullName.Data()
  dataRegion = data_set.reduce(regionCatStr)
  
  nobsRegion = 0.
  
  if dataRegion:
    nobsRegion = dataRegion.sumEntries()
  else:
    print " ERROR : dataset-category dataRegion not found"
    
  if chosenSample:
    regSys['sqrtnobsa'] = 0.
  else:
    regSys['sqrtnobsa'] = TMath.Sqrt(nobsRegion)

  ####

  if chosenSample:
    pdfInRegion  = Util.GetComponent(w,sample,region)
  else:
    rawPdfInRegion = Util.GetRegionPdf(w, region)
    varInRegion =  Util.GetRegionVar(w, region)
    prodList = rawPdfInRegion.pdfList()
    foundRRS = 0
    for idx in range(prodList.getSize()):
      if prodList[idx].InheritsFrom("RooRealSumPdf"):
        rrspdfInt =  prodList[idx].createIntegral(RooArgSet(varInRegion));
        pdfInRegion = rrspdfInt
        foundRRS += 1
    if foundRRS >1 or foundRRS==0:
      print " \n\n WARNING: ", pdf.GetName(), " has ", foundRRS, " instances of RooRealSumPdf"
      print pdf.GetName(), " component list:", prodList.Print("v")
    
  if not pdfInRegion:
    if chosenSample:
      print " \n Warning, could not find pdf in region = ",region, " for sample = ",sample
    else:
      print " \n Warning, could not find pdf in region = ",region

  nFittedInRegion = pdfInRegion.getVal()
  regSys['sqrtnfitted'] = TMath.Sqrt(nFittedInRegion)
  regSys['nfitted'] = nFittedInRegion

  pdfFittedErrInRegion = Util.GetPropagatedError(pdfInRegion, result, doAsym) 
  regSys['totsyserr'] = pdfFittedErrInRegion


  # calculate error per parameter on  fitresult
  fpf = result.floatParsFinal() 
  
  # set all floating parameters constant
  for idx in range(fpf.getSize()):
    parname = fpf[idx].GetName()
    par = w.var(parname)
    par.setConstant()
    
  for idx in range(fpf.getSize()):
    parname = fpf[idx].GetName()
    par = w.var(parname)
    par.setConstant(False)
    sysError  = Util.GetPropagatedError(pdfInRegion, result, doAsym)
    if namemap.has_key(parname): ## add this if I want description
      parname = namemap[parname] ## add this if I want description
    regSys['syserr_'+parname] =  sysError
    par.setConstant() 

  

  return regSys





##################################
##################################
##################################


def latexfitresults_method2(filename,resultname='RooExpandedFitResult_afterFit', region='3jL', sample='', fitregions = 'WR,TR,S3,S4,SR3jT,SR4jT', dataname='obsData', doAsym=False):

#  namemap = {}
#  namemap = getnamemap()

 ############################################
   
  w = Util.GetWorkspaceFromFile(filename,'w')
  if w==None:
    print "ERROR : Cannot open workspace : "
    sys.exit(1) 

  result = w.obj(resultname)
  if result==None:
    print "ERROR : Cannot open fit result : ", resultname
    sys.exit(1)

  resultlistOrig = result.floatParsFinal()
    
  snapshot =  'snapshot_paramsVals_' + resultname
  w.loadSnapshot(snapshot)

  data_set = w.data(dataname)
  if data_set==None:
    print "ERROR : Cannot open dataset : ", "data_set"
    sys.exit(1)
      
  regionCat = w.obj("channelCat")
  data_set.table(regionCat).Print("v");

  regionFullName = Util.GetFullRegionName(regionCat, region)

  fitRegionsList = fitregions.split(",")
  fitRegionsFullName = ""
  for reg in fitRegionsList:
    regFullName = Util.GetFullRegionName(regionCat, reg)
    if fitRegionsFullName == "":
      fitRegionsFullName = regFullName.Data()
    else:
      fitRegionsFullName = fitRegionsFullName + "," + regFullName.Data()

  chosenSample = False
  if sample is not '':
    chosenSample = True

  #####################################################

  regSys = {}

  regionCatStr = 'channelCat==channelCat::' + regionFullName.Data()
  dataRegion = data_set.reduce(regionCatStr)
  nobsRegion = 0.
  
  if dataRegion:
    nobsRegion = dataRegion.sumEntries()
  else:
    print " ERROR : dataset-category", regionCatStr, " not found"
    
  if chosenSample:
    regSys['sqrtnobsa'] = 0.
  else:
    regSys['sqrtnobsa'] = TMath.Sqrt(nobsRegion)

  ####

  if chosenSample:
    pdfInRegion  = Util.GetComponent(w,sample,region)
  else:
    rawPdfInRegion = Util.GetRegionPdf(w, region)
    varInRegion =  Util.GetRegionVar(w, region)
    prodList = rawPdfInRegion.pdfList()
    foundRRS = 0
    for idx in range(prodList.getSize()):
      if prodList[idx].InheritsFrom("RooRealSumPdf"):
        rrspdfInt =  prodList[idx].createIntegral(RooArgSet(varInRegion));
        pdfInRegion = rrspdfInt
        foundRRS += 1
    if foundRRS >1 or foundRRS==0:
      print " \n\n WARNING: ", pdf.GetName(), " has ", foundRRS, " instances of RooRealSumPdf"
      print pdf.GetName(), " component list:", prodList.Print("v")
    
  if not pdfInRegion:
    if chosenSample:
      print " \n Warning, could not find pdf in region = ",region, " for sample = ",sample
    else:
      print " \n Warning, could not find pdf in region = ",region

  nFittedInRegion = pdfInRegion.getVal()
  regSys['sqrtnfitted'] = TMath.Sqrt(nFittedInRegion)

  pdfFittedErrInRegion = Util.GetPropagatedError(pdfInRegion, result, doAsym) 
  regSys['totsyserr'] = pdfFittedErrInRegion
  
  # redo the fit for every parameter being fixed
  lumiConst = True
  fpf = result.floatParsFinal()
  
  # redo the fit for every parameter being fixed
  for idx in range(fpf.getSize()):
    
    parname = fpf[idx].GetName()
    print "\n Method-2: redoing fit with fixed parameter ", parname

    # the parameter that is fixed, needs to have the value of the default fit
    w.loadSnapshot(snapshot)
    par = w.var(parname)

    #     # before redoing the fit, set the values of parameters to initial snapshot, otherwise MIGRAD cannot find improvement
    #     w.loadSnapshot('snapshot_paramsVals_initial')
    #     par.setVal(parDefVal)
    par.setConstant(True)
    suffix = parname + "Fixed"
    result_1parfixed = Util.FitPdf(w, fitRegionsFullName, lumiConst, data_set, suffix)

    expResultAfter_1parfixed = RooExpandedFitResult(result_1parfixed, resultlistOrig)

    nFittedInRegion_1parfixed = pdfInRegion.getVal()
    pdfFittedErrInRegion_1parfixed = Util.GetPropagatedError(pdfInRegion, expResultAfter_1parfixed, doAsym) #  result_1parfixed)

    if pdfFittedErrInRegion_1parfixed > pdfFittedErrInRegion:
      print "\n\n  WARNING  parameter ", parname," gives a larger error when set constant. Do you expect this?"
      print "  WARNING          pdfFittedErrInRegion = ", pdfFittedErrInRegion, "    pdfFittedErrInRegion_1parfixed = ", pdfFittedErrInRegion_1parfixed

    systError  =  TMath.Sqrt(abs(pdfFittedErrInRegion*pdfFittedErrInRegion - pdfFittedErrInRegion_1parfixed*pdfFittedErrInRegion_1parfixed))
    par.setConstant(False)

    if result_1parfixed.status()==0 and result_1parfixed.covQual()==3:   #and result_1parfixed.numStatusHistory()==2 and  result_1parfixed.statusCodeHistory(0)==0 and  result_1parfixed.statusCodeHistory(1) ==0:
      systError = systError
    else:
      systError = 0.0
      print "        WARNING :   for parameter ",parname," fixed the fit does not converge, as status=",result_1parfixed.status(), "(converged=0),  and covariance matrix quality=", result_1parfixed.covQual(), " (full accurate==3)"
      print "        WARNING: setting systError = 0 for parameter ",parname

      #if namemap.has_key(parname):
      #  parname = namemap[parname]
    regSys['syserr_'+parname] =  systError

  return regSys




##################################
##################################
##################################

# MAIN

if __name__ == "__main__":
  
  import os, sys
  import getopt
  def usage():
    print "Usage:"
    print "SysTable.py [-c channels] [-w workspace_afterFit] [-o outputFileName] [-o outputFileName] [-s sample] [-m method] [-f fitregions] [-%] [-b]\n"
    print "Minimal set of inputs [-c channels] [-w workspace_afterFit]"
    print "*** Options are: "
    print "-c <channels>: single channel (region) string or comma separated list accepted (OBLIGATORY)"
    print "-w <workspaceFileName>: single name accepted only (OBLIGATORY) ;   if multiple channels/regions given in -c, assumes the workspace file contains all channels/regions"
    print "-s <sample>: single unique sample name or comma separated list accepted (sample systematics will be calculated for every region given)"
    print "-o <outputFileName>: sets the output table file name, name defined by regions if none provided"
    print "-b: shows the error on samples Before the fit (by default After fit is shown)"
    print "-%: also show the individual errors as percentage of the total systematic error (off by default)"
    print "-y: take symmetrized average of minos errors"

    print "\nFor example:"
    print "SysTable.py -w /afs/cern.ch/user/k/koutsman/HistFitterUser/MET_jets_leptons/results/Combined_KFactorFit_5Channel_Validation_combined_BasicMeasurement_model_afterFit.root  -c SR7jTEl_meffInc,SR7jTMu_meffInc"
    print "SysTable.py -w  /afs/cern.ch/user/c/cote/susy0/users/cote/HistFitter5/results/Combined_KFactorFit_5Channel_bkgonly_combined_BasicMeasurement_model_afterFit.root  -c SR7jTEl_meffInc,SR7jTMu_meffInc -o SystematicsMultiJetsSR.tex"
    print "SysTable.py -w  /afs/cern.ch/user/k/koutsman/HistFitterUser/MET_jets_leptons/results/Combined_KFactorFit_5Channel_Validation_combined_BasicMeasurement_model_afterFit.root  -c SR7jTEl,SR7jTMu -m 2 -f WREl,WRMu,TREl,TRMu"
    print "SysTable.py -w  /afs/cern.ch/user/k/koutsman/HistFitterUser/MET_jets_leptons/results/Combined_KFactorFit_5Channel_Validation_combined_BasicMeasurement_model_afterFit.root  -c SR7jTEl,SR7jTMu -s Top,WZ"
    print "SysTable.py -w ~/Combined_KFactorFit_5Channel_Validation_combined_BasicMeasurement_model_afterFit.root -c SR7jTEl -m 2 -f TRee_nJet,TRem_nJet,TRmm_nJet,TREl_nJet,TRMu_nJet,ZRee_nJet,ZRmm_nJet,WREl_nJet,WRMu_nJet"

    print "\n  Method-1: set all parameters constant, except for the one you're interested in, calculate the error propagated due to that parameter"
    print "  Method-2: set the parameter you're interested in constant, redo the fit with all other parameters floating, calculate the quadratic difference between default fit and your new model with parameter fixed"
    sys.exit(0)        

  wsFileName='/results/MyOneLeptonKtScaleFit_HardLepR17_BkgOnlyKt_combined_NormalMeasurement_model_afterFit.root'
  try:
    opts, args = getopt.getopt(sys.argv[1:], "o:c:w:m:f:s:%by")
  except:
    usage()
  if len(opts)<2:
    usage()

  outputFileName="default"
  method="1"
  showAfterFitError=True
  showPercent=False
  doAsym=True
  sampleStr=''
  for opt,arg in opts:
    if opt == '-c':
      chanStr=arg.replace(",","_")
      chanList=arg.split(",")
    elif opt == '-w':
      wsFileName=arg
    elif opt == '-o':
      outputFileName=arg
    elif opt == '-m':
      if arg == "2" or arg == "1":
        method = arg
      else:
        print "Warning, only methods 1 or 2 are possible. You set method (-m) = ", arg
        sys.exit(0) 
    elif opt == '-f':
      fitRegionsStr=arg
      fitRegionsList=arg.split(",")
    elif opt == '-s':
      sampleStr=arg.replace(",","_") + "_"
      sampleList=arg.split(",")
    elif opt == '-b':
      showAfterFitError=False
    elif opt == '-%':
      showPercent=True
    elif opt == '-y':
      doAsym=True
     
  if outputFileName=="default":
    outputFileName=sampleStr+chanStr+'_SysTable.tex'
    pass

  try:
    fitRegionsList
    if fitRegionsList and not method=="2":
      print "Warning, you set fitRegions (-f) = ", fitRegionsStr, " but not method 2 (-m 2). Fitregions can only be set together with method 2"
      sys.exit(0)
  except NameError:
    pass

  if method=="2":
    try:
      fitRegionsList
    except NameError:
      print "Warning, you did not set fitRegions (-f), but set method 2 (-m 2). Fitregions must be specified when running method 2"
      sys.exit(0)

  chosenSample = False
  try:
    sampleList
    chosenSample=True
  except NameError:
    pass

 
  resultName = 'RooExpandedFitResult_afterFit'
  if not showAfterFitError:
    resultName =  'RooExpandedFitResult_beforeFit'

  skiplist = ['sqrtnobsa', 'totbkgsysa', 'poisqcderr','sqrtnfitted','totsyserr','nfitted']

  chanSys = {}
  origChanList = list(chanList)
  chanList = []

  for chan in origChanList:

    # this 'try' is a hack, I don't care about this software enough to
    # make it good. DG
    try:
      if not chosenSample:
          if method == "2":
              regSys = latexfitresults_method2(wsFileName,resultName,chan,'',fitRegionsStr,'obsData',doAsym)
          else:
              regSys = latexfitresults(wsFileName,chan,'',resultName,'obsData',doAsym)
          chanSys[chan] = regSys
          chanList.append(chan)
      else:
        for sample in sampleList:
          if method == "2":
            regSys = latexfitresults_method2(wsFileName,resultName,chan,sample,fitRegionsStr,'obsData',doAsym)
          else:
            regSys = latexfitresults(wsFileName,chan,sample,resultName,'obsData',doAsym)
          chanSys[chan+"_"+sample] = regSys
          chanList.append(chan+"_"+sample)
    except AttributeError:
      pass

  line_chanSysTight = tablefragment(chanSys,'Signal',chanList,skiplist,chanStr,showPercent)
  
  f = open(outputFileName, 'w')
  f.write( line_chanSysTight )
  f.close()
  print "\nwrote results in file: %s"%(outputFileName)

