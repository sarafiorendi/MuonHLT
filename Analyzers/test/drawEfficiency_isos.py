import numpy
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i"  , "--input"     , dest = "input"     ,  help = "input file"      , default = ''                        )
parser.add_argument("-c"  , "--compare"   , dest = "compfile"  ,  help = "file to compare" , default = ''                        )
parser.add_argument("-d"  , "--diff"      , dest = "diff"      ,  help = "plot differences", default = False, action='store_true')
parser.add_argument("-m"  , "--mc"        , dest = "mc"        ,  help = "comparison is mc", default = False, action='store_true')


options = parser.parse_args()
if not options.input:   
  parser.error('Input filename not given')

import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis, TPad, TGraphErrors
from   array import array
from   math  import sqrt, isnan

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

# file1  = TFile.Open('efficiency_onZMuMuSkim.root'      , 'r')
file = TFile.Open(options.input  , 'r')
if options.compfile:
  file_comp = TFile.Open(options.compfile  , 'r')
# file2  = TFile.Open('efficiency_onDYJetsToLL_MC.root'  , 'r')

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


pt_bins  = [ 20, 24, 27, 30, 33, 36, 39, 42, 45, 49, 53, 57 ,62 ,67 ,72 ,80 , 90, 110, 150] 
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
    pEff1.SetMarkerStyle(22)
    pEff1.SetMarkerSize(0.7)

    pEff1.SetTitle(";" + var[2] + ";" + var[3])
    
#     pEff1.Draw("AP")
#     
#     gPad.Update()
# 
#     pEff1.GetPaintedGraph().GetXaxis().SetLimits(var[5][0], var[5][1])
#     pEff1.GetPaintedGraph().GetHistogram().SetMinimum(var[6][0])             
#     pEff1.GetPaintedGraph().GetHistogram().SetMaximum(var[6][1])        
#   
#     pEff1.GetPaintedGraph().GetXaxis().SetLabelSize(0.04)
#     pEff1.GetPaintedGraph().GetYaxis().SetLabelSize(0.04)   
#     pEff1.GetPaintedGraph().GetXaxis().SetTitleSize(0.04)
#     pEff1.GetPaintedGraph().GetYaxis().SetTitleSize(0.04)
#     pEff1.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
#     pEff1.GetPaintedGraph().GetXaxis().SetTitleOffset(1.1)
#     
    return pEff1




    


atitle = 'offline #Delta#beta corrected rel iso'
ytitle = 'isolation efficiency'

variables = [
#  numerator          # denominator           # x axis title            # y title   # rebin  # x range     # y range  # pdf name                    # legend position                
 ('muonPt_num'       , 'muonPt_den'        , 'muon p_{T} [GeV]'        , ytitle,  pt_bins  , (  20  , 150), (0.8, 1),  'efficiency_muonPt'       ,  (0.6 , 0.85, 0.37, 0.52), (0.901, 1.1 )), 
 ('muonEta_num'      , 'muonEta_den'       , 'muon #eta '              , ytitle,  eta_bins , ( -2.4 , 2.4), (0.8, 1),  'efficiency_muonEta'      ,  (0.6 , 0.85, 0.37, 0.52), (0.901, 1.1 )),
 ('muonIso_num'      , 'muonIso_den'       ,  atitle + ', #DeltaR=0.3' , ytitle,  iso_bins , (   0  , 0.6), (0. , 1),  'efficiency_muonIso03'    ,  (0.69, 0.88, 0.70, 0.85), (0.   , 2   )),
 ('muonIso04_num'    , 'muonIso04_den'     ,  atitle + ', #DeltaR=0.4' , ytitle,  iso_bins , (   0  , 0.6), (0. , 1),  'efficiency_muonIso04'    ,  (0.69, 0.88, 0.70, 0.85), (0.   , 2   )),
 ('muonPt_barrel_num', 'muonPt_barrel_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins  , (  20  , 150), (0.8, 1),  'efficiency_muonPt_barrel',  (0.6 , 0.85, 0.37, 0.52), (0.901, 1.1 )),
 ('muonPt_endcap_num', 'muonPt_endcap_den' , 'muon p_{T} [GeV]'        , ytitle,  pt_bins  , (  20  , 150), (0.8, 1),  'efficiency_muonPt_endcap',  (0.6 , 0.85, 0.37, 0.52), (0.901, 1.1 )),
 ('nvtx_num'         , 'nvtx_den'          , '# offline PV'            , ytitle,         1 , (   3  ,  31), (0.7, 1),  'efficiency_NPV'          ,  (0.5 , 0.78, 0.37, 0.52), (0.6  , 1.4 )),
 ('nvtx_barrel_num'  , 'nvtx_barrel_den'   , '# offline PV'            , ytitle,         1 , (   3  ,  31), (0.7, 1),  'efficiency_NPV_barrel'   ,  (0.5 , 0.78, 0.37, 0.52), (0.6  , 1.4 )),
 ('nvtx_endcap_num'  , 'nvtx_endcap_den'   , '# offline PV'            , ytitle,         1 , (   3  ,  31), (0.7, 1),  'efficiency_NPV_endcap'   ,  (0.5 , 0.78, 0.37, 0.52), (0.6  , 1.4 )),
 ('muonPhi_num'      , 'muonPhi_den'       , 'muon #phi '              , ytitle,         1 , ( -3.2 , 3.2), (0.8, 1),  'efficiency_muonPhi'      ,  (0.5 , 0.78, 0.37, 0.52), (0.9  , 1.1 )),
] 


for var in variables:
  pEff1 = doHisto(file , var, ROOT.kBlue)

#   c.cd()
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


  if options.compfile:
    if not options.mc:
      pEff2 = doHisto(file_comp , var, ROOT.kRed)
    else:
#       pEff2 = file_comp.Get(var[0].split('_')[0]  )
      pEff2 = file_comp.Get(var[0].replace('_num',''))
      pEff2.SetLineColor  (ROOT.kRed)
      pEff2.SetMarkerColor(ROOT.kRed)
      pEff2.SetMarkerStyle(22)
      pEff2.SetMarkerSize(0.7)

    pEff2.Draw('P same')
    c.Update()
    c.Modified()
#     pEff3
#     pEff3.Draw('P same')
    l = TLegend(var[8][0], var[8][2], var[8][1], var[8][3])
    l.SetBorderSize(0)
    l.SetTextSize(0.037)
    l.AddEntry(pEff1 , "Data"   , "pel")
    l.AddEntry(pEff2 , "MC"     , "pel")
    if var[0].find('barrel') > 0:  
      l.SetHeader('barrel')
    if var[0].find('endcap') > 0:  
      l.SetHeader('endcap')
    l.Draw()
    
    # draw bottom pad with ratio h1/h2
    if options.diff:
      c.SetCanvasSize(600,800)
#       c.SetCanvasSize(600,900)
      stackPad.SetBottomMargin(0.2)
      ratioPad.cd()
      ratioPad.SetGridy(True)
      ratioPad.SetBottomMargin(0.2)
      h_base_clone = pEff1.GetPaintedGraph()
      h_comp_clone = pEff2.GetPaintedGraph()
      h_for_errors = pEff1.GetPassedHistogram()
      N            = h_base_clone.GetN()
      bin_base     = numpy.frombuffer(h_base_clone.GetY(),  count = N)
      bin_comp     = numpy.frombuffer(h_comp_clone.GetY(),  count = N)
      bin_x        = numpy.frombuffer(h_comp_clone.GetX(),  count = N)

      bin_ey_l_b   = numpy.frombuffer(h_base_clone.GetEYlow(),  count = N)
      bin_ey_h_b   = numpy.frombuffer(h_base_clone.GetEYhigh(), count = N)
      bin_ey_max_b = numpy.maximum(bin_ey_l_b, bin_ey_h_b) 

      bin_ey_l_c   = numpy.frombuffer(h_comp_clone.GetEYlow(),  count = N)
      bin_ey_h_c   = numpy.frombuffer(h_comp_clone.GetEYhigh(), count = N)
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
#       bin_eratio_tmp = [( a/b*sqrt((ea*ea)/(a*a) + (eb*eb)/(b*b)) ) for a,b,ea,eb in zip(bin_base,bin_comp,bin_ey_max_b,bin_ey_max_c) if a!=0 and b!=0 and eb!=0]
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
#       import pdb; pdb.set_trace()
      ratioGraph = TGraphErrors(N, bin_x, vec_ratio, vec_ex, vec_eratio)
      ratioGraph.SetTitle('')
      ratioGraph.GetYaxis().SetTitle('data/MC')
      ratioGraph.GetXaxis().SetTitle(pEff1.GetPaintedGraph().GetXaxis().GetTitle())
      ratioGraph.Draw("AP")
      ratioGraph.GetYaxis().SetRangeUser( var[9][0], var[9][1])
      ratioGraph.GetXaxis().SetRangeUser( var[5][0], var[5][1])#var[4][0], var[4][1])
      ratioGraph.GetYaxis().SetLabelSize(0.09)
      ratioGraph.GetXaxis().SetLabelSize(0.1 )
      ratioGraph.GetYaxis().SetTitleSize(0.08)
      ratioGraph.GetXaxis().SetTitleSize(0.09)  
      ratioGraph.GetYaxis().SetTitleOffset(0.8)
      ratioGraph.GetYaxis().SetNdivisions(505)
#       h_ratio.Draw('e')


  gPad.SetGridx(True)
  gPad.SetGridy(True)
  c.SaveAs(var[7] + "_ZMuMuSkim_run251643_comparison_hltisoTag_offlineIsoTag_match.pdf")







