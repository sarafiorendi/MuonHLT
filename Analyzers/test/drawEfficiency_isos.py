import numpy
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i"  , "--input"     , dest = "input"     ,  help = "input file"       , default = ''                        )
parser.add_argument("-c"  , "--compare"   , dest = "compfile"  ,  help = "file to compare"  , default = ''                        )
parser.add_argument("-d"  , "--diff"      , dest = "diff"      ,  help = "plot differences" , default = False, action='store_true')
parser.add_argument("-m"  , "--mc"        , dest = "mc"        ,  help = "comparison is mc" , default = False, action='store_true')
parser.add_argument("-l"  , "--leg"       , dest = "leg"       ,  help = "legend labels"    , default = ''                        )
parser.add_argument("-j"  , "--input2"    , dest = "input2"    ,  help = "input file"       , default = ''                        )




options = parser.parse_args()
if not options.input:   
  parser.error('Input filename not given')

import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis, TPad, TGraphErrors
from   array import array
from   math  import sqrt, isnan
from   copy  import deepcopy 

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

# import pdb; pdb.set_trace()

namefiles = options.input.split(',')
nfiles   = len(namefiles)
files    = []

print 'number of input files is ' + str(nfiles)

for i in range(0, nfiles):
  print 'opening file ' + str(i) + ': ' + namefiles[i]
  files.append(TFile.Open(namefiles[i]   , 'r') )


c = TCanvas('', '', 600,600)
if options.diff:
  stackPad = ROOT.TPad('stackPad', 'stackPad', 0.,  .25, 1., 1.  , 0, 0)  
  ratioPad = ROOT.TPad('ratioPad', 'ratioPad', 0., 0. , 1.,  .3, 0, 0)  
else:
  stackPad = ROOT.TPad('stackPad', 'stackPad', 0,  0. , 1., 1.  , 0, 0)  
  ratioPad = ROOT.TPad('ratioPad', 'ratioPad', 0., 0. , 0., 0.  , 0, 0)  
c.cd()
stackPad.Draw()
ratioPad.Draw()

pt_bins  = [  0, 15, 18, 20, 22, 25, 30, 40, 50, 60, 80, 120] 
eta_bins = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]

colorlist = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kAzure+1, ROOT.kViolet]


def doHisto(file, var, thecolor, i):

  if i==0:
    pEff1  = file.Get(var[0] + '2016'  )
  else: 
    pEff1  = file.Get(var[0]  )
#   pEff1  = file.Get(var[0]  )
    
  pEff1.SetLineColor  (thecolor)
  pEff1.SetMarkerColor(thecolor)
  pEff1.SetMarkerStyle(8  )
  pEff1.SetMarkerSize(0.8)

  pEff1.SetTitle(";" + var[1] + ";" + var[2])
  return pEff1





# ytitle = 'TkMu/L1 efficiency'
ytitle = 'trk isolation efficiency'

variables = [
#  numerator          # x axis title            # y title   # rebin    # x range      # y range      # pdf name                     # legend position         #y range ratio          
#  ('muonPt_barrel'    , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0. , 1.01),  'efficiency_muonPt_barrel_L2',  (0.3 , 0.75, 0.18, 0.3), (0.901, 1.05 )), 
#  ('muonPt_endcap'    , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0. , 1.01),  'efficiency_muonPt_endcap_L2',  (0.3 , 0.75, 0.18, 0.3), (0.901, 1.05 )), 
#  ('muonPt_barrel'    , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0.8, 1.01),  'efficiency_muonPt_barrel_L2_zoom',  (0.34 , 0.8, 0.18, 0.35), (0.901, 1.05 )), 
#  ('muonPt_endcap'    , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0.8, 1.01),  'efficiency_muonPt_endcap_L2_zoom',  (0.34 , 0.8, 0.18, 0.35), (0.901, 1.05 )), 
#  ('muonPt_BMTF'      , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0.8, 1.01),  'efficiency_muonPt_BMTF_L2'  ,  (0.3 , 0.75, 0.18, 0.35), (0.901, 1.05 )), 
#  ('muonPt_OMTF'      , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0.3, 1.01),  'efficiency_muonPt_OMTF_L2'  ,  (0.3 , 0.75, 0.18, 0.35), (0.901, 1.05 )), 
#  ('muonPt_EMTF'      , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0   , 150), (0.2, 1.01),  'efficiency_muonPt_EMTF_L2'  ,  (0.3 , 0.75, 0.18, 0.35), (0.901, 1.05 )), 
 ('muonPt'           , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  0    , 150 ), (0.85 , 1.01),  'efficiency_muonPt'   ,  (0.3 , 0.75, 0.18, 0.35), (0.901, 1.05 )), 
 ('muonEta'          , 'muon #eta '              , ytitle,   eta_bins , ( -2.4  , 2.4 ), (0.95 , 1.0),  'efficiency_muonEta'  ,  (0.3 , 0.6, 0.18, 0.32),  (0.9  , 1.05  )),
 ('muonPhi'          , 'muon #phi '              , ytitle,   1        , ( -3.14 , 3.14), (0.95 , 1.0),  'efficiency_muonPhi'  ,  (0.3 , 0.6, 0.18, 0.32),  (0.9  , 1.05  )),
] 



for var in variables:
  
  l = TLegend(var[7][0], var[7][2], var[7][1], var[7][3])
  l.SetBorderSize(0)
  l.SetTextSize(0.028)

  for i, ifile in enumerate(files):
    pEff1 = doHisto(ifile , var, colorlist[i], i)

    stackPad.cd()
    if (i == 0):
      pEff1.Draw('AP')
      c.Update()
      c.Modified()
      pEff1.GetPaintedGraph().GetXaxis().SetLimits(var[4][0], var[4][1])
      pEff1.GetPaintedGraph().GetHistogram().SetMinimum(var[5][0])             
      pEff1.GetPaintedGraph().GetHistogram().SetMaximum(var[5][1])        
  
      pEff1.GetPaintedGraph().GetXaxis().SetLabelSize(0.04)
      pEff1.GetPaintedGraph().GetYaxis().SetLabelSize(0.04)   
      pEff1.GetPaintedGraph().GetXaxis().SetTitleSize(0.04)
      pEff1.GetPaintedGraph().GetYaxis().SetTitleSize(0.04)
      pEff1.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
      pEff1.GetPaintedGraph().GetXaxis().SetTitleOffset(1.2)

    else:
      pEff1.Draw('P same')

    c.Update()
    c.Modified()

    if options.leg:
#       labels = []
#       labels.append( options.leg.split(',')[i] )
      l.AddEntry(pEff1 , options.leg.split(',')[i]  , "pel")


#   if var[0].find('barrel') > 0:  
#     l.SetHeader('barrel')
#   if var[0].find('endcap') > 0:  
#     l.SetHeader('endcap')
#   if var[0].find('eta0') > 0:  
#     l.SetHeader('0 < |#eta| < 0.9')
#   if var[0].find('eta1') > 0:  
#     l.SetHeader('0.9 < |#eta| < 1.2')
#   if var[0].find('eta2') > 0:  
#     l.SetHeader('1.2 < |#eta| < 2.4')

  l.Draw()
    


  gPad.SetGridx(True)
  gPad.SetGridy(True)
  c.SaveAs("" +  var[6] + "_TkIso2017.pdf")



