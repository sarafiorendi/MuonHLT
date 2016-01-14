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

file  = TFile.Open(options.input   , 'r')
file2 = TFile.Open(options.input2  , 'r')
if options.compfile:
  file_comp   = TFile.Open(options.compfile.split(',')[0]   , 'r')
  file_comp2  = TFile.Open(options.compfile.split(',')[1]   , 'r')



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


pt_bins       = [  18, 20, 22, 25, 30, 40, 50, 60, 90, 150] 
# pt_bins_less  = [ 15, 18, 24, 30, 40, 50, 60, 90, 150]
# pt_bins  = [ 15, 18, 20, 24, 27, 30, 35, 40, 45, 50, 60 , 70 , 90, 150] 
# pt_bins  = [ 18, 20, 24, 27, 30, 33, 36, 39, 42, 45, 49, 53, 57 ,62 ,67 ,72 ,80 , 90, 110, 150] 
eta_bins = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
iso_bins = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1]
 


def doHisto(file, var, thecolor):
  h_num  = file.Get(var[0]  )
  h_den  = file.Get(var[1]  )
  
  if isinstance(var[4], list): 
    thebin = array ('d',var[4])
    h_num_r = h_num.Rebin( len(var[4])-1, 'h_num_r', thebin)
    h_den_r = h_den.Rebin( len(var[4])-1, 'h_den_r', thebin)
  else:   
    h_num_r = h_num.Rebin(var[4], 'h_num_r')
    h_den_r = h_den.Rebin(var[4], 'h_den_r')

  h_lowbin = h_num_r.GetXaxis().GetBinLowEdge(1)  
  h_upbin  = h_num_r.GetXaxis().GetBinUpEdge(h_num_r.GetNbinsX())  
  pEff1    = TEfficiency('eff1',"my efficiency",h_num_r.GetNbinsX(), h_lowbin, h_upbin)

  if TEfficiency.CheckConsistency( h_num_r, h_den_r ):
    pEff1 = TEfficiency(h_num_r,h_den_r)
    pEff1.SetLineColor  (thecolor)
    pEff1.SetMarkerColor(thecolor)
    pEff1.SetMarkerStyle(20)
    pEff1.SetMarkerSize(0.8)

    pEff1.SetTitle(";" + var[2] + ";" + var[3])
    return pEff1




def doRatio(h_data_clone, h_mc_clone, h_for_errors):

  N            = h_data_clone.GetN()  
  bin_base     = numpy.frombuffer(h_data_clone.GetY(),  count = N)
  bin_comp     = numpy.frombuffer(h_mc_clone.GetY()  ,  count = N)
  bin_x        = numpy.frombuffer(h_data_clone.GetX(),  count = N)

  bin_ey_l_b   = numpy.frombuffer(h_data_clone.GetEYlow(),  count = N)
  bin_ey_h_b   = numpy.frombuffer(h_data_clone.GetEYhigh(), count = N)
  bin_ey_max_b = numpy.maximum(bin_ey_l_b, bin_ey_h_b) 

  bin_ey_l_c   = numpy.frombuffer(h_mc_clone.GetEYlow() , count = N)
  bin_ey_h_c   = numpy.frombuffer(h_mc_clone.GetEYhigh(), count = N)
  bin_ey_max_c = numpy.maximum(bin_ey_l_c, bin_ey_h_c) 
  
    
  bin_eratio_tmp = []
  for i,a in enumerate(bin_base):
    b  = bin_comp[i]
    ea = bin_ey_max_b[i]
    eb = bin_ey_max_c[i]
    if a!=0 and b!=0 and eb!=0:
      bin_eratio_tmp.append(b/a*sqrt((ea*ea)/(a*a) + (eb*eb)/(b*b)) ) 
    else:
      bin_eratio_tmp.append(0)   
    for index, item in enumerate(bin_eratio_tmp):
      if isnan(item):
        bin_eratio_tmp[index] = 0

    bin_eratio  =  numpy.around(bin_eratio_tmp, decimals=5)
    bin_ex      = []
    bin_ratio   = []
    
    for i in range(len(bin_base)):
      if bin_base[i] != 0 and bin_comp[i] != 0: 
        bin_ratio .append ( bin_base[i] / max(0.0000001, bin_comp[i] ))
      else:
        bin_ratio .append(0)

      bin_ex    .append(h_for_errors.GetBinWidth(i+1)/2)
    
    vec_ex     = numpy.asarray(bin_ex    )  
    vec_ratio  = numpy.asarray(bin_ratio )  
    vec_eratio = numpy.asarray(bin_eratio)  

  ratioGraph = TGraphErrors(N, bin_x, vec_ratio, vec_ex, vec_eratio)
  return ratioGraph
    



ytitle = 'efficiency'

variables = [
#  numerator          # denominator          # x axis title            # y title   # rebin    # x range      # y range      # pdf name                  # legend position         #y range ratio          
 ('muonPt_num'       , 'muonPt_den'        , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.6, 1.01),  'efficiency_muonPt'      ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.05 )), 
 ('muonEta_num'      , 'muonEta_den'       , 'muon #eta '              , ytitle,   eta_bins , ( -2.4 , 2.4), (0.6, 1.01),  'efficiency_muonEta'     ,  (0.5 , 0.85, 0.25, 0.45), (0.9  , 1.05  )),
 ('muonPt_eta0_num'  , 'muonPt_eta0_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.6, 1.01),  'efficiency_muonPt_eta0' ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )),
 ('muonPt_eta1_num'  , 'muonPt_eta1_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.6, 1.01),  'efficiency_muonPt_eta1' ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )),
 ('muonPt_eta2_num'  , 'muonPt_eta2_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.6, 1.01),  'efficiency_muonPt_eta2' ,  (0.5 , 0.85, 0.25, 0.45), (0.801, 1.05 )),
 ('nvtx_num'         , 'nvtx_den'          , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.6, 1.01),  'efficiency_NPV'         ,  (0.3 , 0.65, 0.25, 0.45), (0.9  , 1.1  )),
 ('nvtx_eta0_num'    , 'nvtx_eta0_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.6, 1.01),  'efficiency_NPV_eta0'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),
 ('nvtx_eta1_num'    , 'nvtx_eta1_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.6, 1.01),  'efficiency_NPV_eta1'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),
 ('nvtx_eta2_num'    , 'nvtx_eta2_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.6, 1.01),  'efficiency_NPV_eta2'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),  
# these are ranges for fullpath efficiency, at the end of the macro there are the ones for isolation efficiency only

#  ('muonPt_barrel_num', 'muonPt_barrel_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins      , (  15  , 150), (0.6,  1.01),  'efficiency_muonPt_barrel',  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1 )),
#  ('muonPt_endcap_num', 'muonPt_endcap_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins      , (  15  , 150), (0.6,  1.01),  'efficiency_muonPt_endcap',  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1 )),
#  ('nvtx_barrel_num'  , 'nvtx_barrel_den'   , '# offline PV'            , ytitle,         2     , (   2  ,  26), (0.6,  1.01),  'efficiency_NPV_barrel'   ,  (0.3 , 0.65, 0.25, 0.45), (0.6  , 1.4 )),
#  ('nvtx_endcap_num'  , 'nvtx_endcap_den'   , '# offline PV'            , ytitle,         2     , (   2  ,  26), (0.6,  1.01),  'efficiency_NPV_endcap'   ,  (0.3 , 0.65, 0.25, 0.45), (0.6  , 1.4 )),
#  ('muonPhi_num'      , 'muonPhi_den'       , 'muon #phi '              , ytitle,         2     , ( -3.2 , 3.2), (0.6,  1.01),  'efficiency_muonPhi'      ,  (0.5 , 0.85, 0.32, 0.52), (0.9  , 1.1 )),
] 


for var in variables:
  pEff1 = doHisto(file , var, ROOT.kBlack)

  stackPad.cd()
  pEff1.Draw('AP')
  c.Update()
  c.Modified()

  pEff1.GetPaintedGraph().GetXaxis().SetLimits(var[5][0], var[5][1])
  pEff1.GetPaintedGraph().GetHistogram().SetMinimum(var[6][0])             
  pEff1.GetPaintedGraph().GetHistogram().SetMaximum(var[6][1])        

  pEff1.GetPaintedGraph().GetXaxis().SetLabelSize(0.04)
  pEff1.GetPaintedGraph().GetYaxis().SetLabelSize(0.04)   
  pEff1.GetPaintedGraph().GetXaxis().SetTitleSize(0.04)
  pEff1.GetPaintedGraph().GetYaxis().SetTitleSize(0.04)
  pEff1.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
  pEff1.GetPaintedGraph().GetXaxis().SetTitleOffset(1.1)


  if options.input2:
    pEff3 = doHisto(file2 , var, ROOT.kGray+3)
    pEff3.SetMarkerStyle(24)
    pEff3.Draw('P same')
  if options.compfile:
    if not options.mc:
      pEff2 = doHisto(file_comp , var, ROOT.kRed+1)
    else:
      pEff2 = file_comp.Get(var[0].replace('_num',''))
      pEff2.SetLineColor  (ROOT.kRed+1)
      pEff2.SetMarkerColor(ROOT.kRed+1)
      pEff2.SetMarkerStyle(22)
      pEff2.SetMarkerSize(0.8)
      pEff4 = file_comp2.Get(var[0].replace('_num',''))
      pEff4.SetLineColor  (ROOT.kRed)
      pEff4.SetMarkerColor(ROOT.kRed)
      pEff4.SetMarkerStyle(26)
      pEff4.SetMarkerSize(0.8)

    pEff2.Draw('P same')
    pEff4.Draw('P same')
    pEff1.Draw('P same')
    c.Update()
    c.Modified()

    if options.leg:
      l = TLegend(var[8][0], var[8][2], var[8][1], var[8][3])
      l.SetBorderSize(0)
      l.SetTextSize(0.028)
      label1 = options.leg.split(',')[0]
      label2 = options.leg.split(',')[1]
      label3 = options.leg.split(',')[2]
      label4 = options.leg.split(',')[3]
      l.AddEntry(pEff1 , label1  , "pel")
      l.AddEntry(pEff3 , label3  , "pel")
      l.AddEntry(pEff2 , label2  , "pel")
      l.AddEntry(pEff4 , label4  , "pel")
      if var[0].find('barrel') > 0:  
        l.SetHeader('barrel')
      if var[0].find('endcap') > 0:  
        l.SetHeader('endcap')
      if var[0].find('eta0') > 0:  
        l.SetHeader('0 < |#eta| < 0.9')
      if var[0].find('eta1') > 0:  
        l.SetHeader('0.9 < |#eta| < 1.2')
      if var[0].find('eta2') > 0:  
        l.SetHeader('1.2 < |#eta| < 2.4')
      l.Draw()
    
    # draw bottom pad with ratio h1/h2
    if options.diff:
      c.SetCanvasSize(600,800)
      stackPad.SetBottomMargin(0.2)
      ratioPad.cd()
      ratioPad.SetGridy(True)
      ratioPad.SetBottomMargin(0.2)
      h_data_clone = pEff1.GetPaintedGraph()
      h_mc_clone   = pEff2.GetPaintedGraph()
      h_for_errors = pEff1.GetPassedHistogram()

      ratio_isomu = doRatio(h_data_clone, h_mc_clone, h_for_errors)
      ratio_isomu.SetTitle('')
      ratio_isomu.GetYaxis().SetTitle('data/MC')
      ratio_isomu.GetXaxis().SetTitle(pEff1.GetPaintedGraph().GetXaxis().GetTitle())
      ratio_isomu.SetMarkerStyle(21)
      ratio_isomu.SetMarkerSize(0.6)
      ratio_isomu.Draw("AP")
      ratio_isomu.GetYaxis().SetRangeUser( var[9][0], var[9][1])
      ratio_isomu.GetXaxis().SetRangeUser( var[5][0], var[5][1])
      ratio_isomu.GetYaxis().SetLabelSize(0.09)
      ratio_isomu.GetXaxis().SetLabelSize(0.1 )
      ratio_isomu.GetYaxis().SetTitleSize(0.08)
      ratio_isomu.GetXaxis().SetTitleSize(0.09)  
      ratio_isomu.GetYaxis().SetTitleOffset(0.8)
      ratio_isomu.GetYaxis().SetNdivisions(505)

      h_data_clone_old = pEff3.GetPaintedGraph()
      h_mc_clone_old   = pEff4.GetPaintedGraph()
      h_for_errors_old = pEff3.GetPassedHistogram()
      ratio_oldisomu   = doRatio(h_data_clone_old, h_mc_clone_old, h_for_errors_old)
      ratio_oldisomu.SetMarkerStyle(25)
      ratio_oldisomu.SetMarkerSize(0.6)      
      ratio_oldisomu.Draw('P same')

#       l2 = TLegend(var[8][0], 0.6, var[8][1], 0.85)
#       l2.SetBorderSize(0)
#       l2.SetTextSize(0.038)
#       l2.AddEntry(ratio_isomu    , 'HLT_IsoMu18'  , "pel")
#       l2.AddEntry(ratio_oldisomu , 'HLT_OldIsoMu18'  , "pel")
#       l2.Draw()
    
  gPad.SetGridx(True)
  gPad.SetGridy(True)
  c.SaveAs("full" +  var[7] + "_run258158_IsoMu18_offlineIso0p15.pdf")







## for isolation efficiency
# variables = [
# #  numerator          # denominator           # x axis title            # y title   # rebin   # x range     # y range       # pdf name                  # legend position         #y range ratio          
#  ('muonPt_num'       , 'muonPt_den'        , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.75, 1.01),  'efficiency_muonPt'      ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )), 
#  ('muonEta_num'      , 'muonEta_den'       , 'muon #eta '              , ytitle,   eta_bins , ( -2.4 , 2.4), (0.85, 1.01),  'efficiency_muonEta'     ,  (0.5 , 0.85, 0.25, 0.45), (0.95 , 1.05 )),
#  ('muonPt_eta0_num'  , 'muonPt_eta0_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.5,  1.01),  'efficiency_muonPt_eta0' ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )),
#  ('muonPt_eta1_num'  , 'muonPt_eta1_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.7,  1.01),  'efficiency_muonPt_eta1' ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )),
#  ('muonPt_eta2_num'  , 'muonPt_eta2_den'   , 'muon p_{T} [GeV]'        , ytitle,   pt_bins  , (  15  , 150), (0.8,  1.01),  'efficiency_muonPt_eta2' ,  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1  )),
#  ('nvtx_num'         , 'nvtx_den'          , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.8,  1.01),  'efficiency_NPV'         ,  (0.3 , 0.65, 0.25, 0.45), (0.9  , 1.1  )),
#  ('nvtx_eta0_num'    , 'nvtx_eta0_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.75, 1.01),  'efficiency_NPV_eta0'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),
#  ('nvtx_eta1_num'    , 'nvtx_eta1_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.75, 1.01),  'efficiency_NPV_eta1'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),
#  ('nvtx_eta2_num'    , 'nvtx_eta2_den'     , '# offline PV'            , ytitle,       2    , (   2  ,  26), (0.75, 1.01),  'efficiency_NPV_eta2'    ,  (0.3 , 0.65, 0.25, 0.45), (0.8  , 1.2  )),  
# #  ('muonPt_barrel_num', 'muonPt_barrel_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins      , (  15  , 150), (0.6,  1.01),  'efficiency_muonPt_barrel',  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1 )),
# #  ('muonPt_endcap_num', 'muonPt_endcap_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins      , (  15  , 150), (0.6,  1.01),  'efficiency_muonPt_endcap',  (0.5 , 0.85, 0.25, 0.45), (0.901, 1.1 )),
# #  ('nvtx_barrel_num'  , 'nvtx_barrel_den'   , '# offline PV'            , ytitle,         2     , (   2  ,  26), (0.6,  1.01),  'efficiency_NPV_barrel'   ,  (0.3 , 0.65, 0.25, 0.45), (0.6  , 1.4 )),
# #  ('nvtx_endcap_num'  , 'nvtx_endcap_den'   , '# offline PV'            , ytitle,         2     , (   2  ,  26), (0.6,  1.01),  'efficiency_NPV_endcap'   ,  (0.3 , 0.65, 0.25, 0.45), (0.6  , 1.4 )),
# #  ('muonPhi_num'      , 'muonPhi_den'       , 'muon #phi '              , ytitle,         2     , ( -3.2 , 3.2), (0.6,  1.01),  'efficiency_muonPhi'      ,  (0.5 , 0.85, 0.32, 0.52), (0.9  , 1.1 )),
# ] 
