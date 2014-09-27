#!/usr/bin/env python2.7
"""
I've tried to clean this up and understand it, but it's pretty clear
the author couldn't give a flying fuck about writing good code, so not
a huge effort...
"""

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import gROOT
from scharmfit.utils import load_susyfit
load_susyfit()
gROOT.Reset()

from ROOT import Util, TMath, RooStats, RooArgSet
from ROOT import RooFit

from UpperLimitTableTex import tablefragment

import argparse
import os
import sys

## MB instructions
# this is something I wrote, badly, and now I expect you to figure it
# out. Enjoy!

def latexfitresults(filename, poiname='mu_Sig', lumiFB=1.0,
                    nTOYS=3000, asimov=False, wname='combined'):

  workspacename=wname

  w = Util.GetWorkspaceFromFile(filename,workspacename)

  if w==None:
    print "ERROR : Cannot open workspace : ", workspacename
    sys.exit(1)

  if len(poiname)==0:
    print " "
  else:
    modelConfig = w.obj("ModelConfig")
    poi = w.var(poiname)
    if poi==None:
      print "ERROR : Cannot find POI with name: ", poiname, " in workspace from file ", filename
      sys.exit(1)
    modelConfig.SetParametersOfInterest(RooArgSet(poi))
    modelConfig.GetNuisanceParameters().remove(poi)

  ntoys = 3000
  calctype = 0   # toys = 0, asymptotic (asimov) = 2
  npoints = 20

  if nTOYS != 3000 and nTOYS>0:
    ntoys = nTOYS
  if asimov:
    calctype = 2

  hti_result = RooStats.MakeUpperLimitPlot(
    poiname,w,calctype,3,ntoys,True,npoints)
  outFileName = "./htiResult_poi_" + poiname + "_ntoys_" + str(ntoys) + "_calctype_" + str(calctype) + "_npoints_" + str(npoints) + ".root"
  hti_result.SaveAs(outFileName)
  hti_result.Print()

  uL_nobsinSR = hti_result.UpperLimit()
  uL_visXsec = uL_nobsinSR / lumiFB
  # uL_visXsecErrorUp = uL_visXsec - uL_nobsinSR/(lumiFB * (1. + lumiRelUncert))
  # uL_visXsecErrorDown = uL_nobsinSR/(lumiFB * (1. - lumiRelUncert)) - uL_visXsec

  uL_nexpinSR = hti_result.GetExpectedUpperLimit(0)
  uL_nexpinSR_P = hti_result.GetExpectedUpperLimit(1) 
  uL_nexpinSR_M = hti_result.GetExpectedUpperLimit(-1)
  if uL_nexpinSR > uL_nexpinSR_P or uL_nexpinSR < uL_nexpinSR_M:
    print " \n something very strange, either the uL_nexpinSR > uL_nexpinSR_P or uL_nexpinSR < uL_nexpinSR_M"
    print "  uL_nexpinSR = ", uL_nexpinSR , " uL_nexpinSR_P = ", uL_nexpinSR_P, " uL_nexpinSR_M = ", uL_nexpinSR_M
  uL_nexpinSRerrP = hti_result.GetExpectedUpperLimit(1) - uL_nexpinSR
  uL_nexpinSRerrM = uL_nexpinSR - hti_result.GetExpectedUpperLimit(-1)

  # find the CLB values at indexes above and below observed CLs p-value
  CLB_P = 0.
  CLB_M = 0.
  mu_P = 0.
  mu_M = 0.
  index_P = 0
  indexFound = False
  for iresult in range(hti_result.ArraySize()):
    xval = hti_result.GetXValue(iresult) 
    yval = hti_result.GetYValue(iresult)
    if xval>uL_nobsinSR and not indexFound:
      index_P = iresult
      CLB_P = hti_result.CLb(iresult)
      mu_P = xval
      if iresult>0:
        CLB_M = hti_result.CLb(iresult-1)
        mu_M = hti_result.GetXValue(iresult-1)
        indexFound = True
 #       print " \n   found the CLB values to interpolate"
 #       print " CLB_M =", CLB_M, " CLB_P =", CLB_P, "  mu_P = ", mu_P, " mu_M = ", mu_M

  # interpolate the value of CLB to be exactly above upperlimit p-val
  alpha_CLB = (CLB_P - CLB_M) / (mu_P - mu_M)
  beta_CLB = CLB_P - alpha_CLB*mu_P
  # CLB is taken as the point on the CLB curve for the same poi value,
  # as the observed upperlimit
  CLB = alpha_CLB * uL_nobsinSR + beta_CLB
  #print " CLB = " , CLB

  print "\n\n\n\n  ***---  now doing p-value calculation ---*** \n\n\n\n"
  pval = RooStats.get_Presult(w,False,1000,2)

  print "p-value is: ", pval
  
  ulList = [uL_visXsec, uL_nobsinSR, uL_nexpinSR, uL_nexpinSRerrP, uL_nexpinSRerrM, CLB, pval ]

  return ulList





##################################
##################################
##################################

#### Main function calls start here ....

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('workspaces', nargs='+')
  parser.add_argument('-o','--output-file-name', default='default')
  parser.add_argument('-n','--n-toys', type=int, default=0)
  args = parser.parse_args()

  outputFileName = "default"
  useAsimovSet = (args.n_toys == 0)
  nTOYS = args.n_toys
  wsFileNameList = args.workspaces
  outputFileName = args.output_file_name

  if outputFileName == "default":
    outputFileName = "UpperLimitTable.tex"
    if useAsimovSet:
      outputFileName = "UpperLimitTable_asimov.tex"
    elif nTOYS >0:
      outputFileName = "UpperLimitTable_nToys"+str(nTOYS)+".tex"

  for tmp_wsFileName in wsFileNameList:
    if not os.path.isfile(tmp_wsFileName):
      print " \n\n\n Warning: workspace file ", tmp_wsFileName, " does not exist, remove this channel from command or make this workspace file available"
      sys.exit(0)

  upLim = {}
  for idx, wsFileNameChan in enumerate(wsFileNameList):
    # calculate upper limit
    ulMapChan = latexfitresults(wsFileNameChan, 'mu_Sig', 20.3, nTOYS, useAsimovSet )#, chan)
    upLim[wsFileNameChan.split('/')[-2]] = ulMapChan

  tablename = 'upperlimit.'
  line_upLim = tablefragment(upLim,tablename)
  f = open(outputFileName,'w')
  f.write( line_upLim )
  f.close()
  print "\nResult written in:"
  print outputFileName

