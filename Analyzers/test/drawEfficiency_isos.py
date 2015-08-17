import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis
from   array import array

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

file1  = TFile.Open('provaOut.root'  , 'r')
# file20_tight = TFile.Open('efficiencyReco_PF_isoMu20_TightMu_offlineRelIso_addNPrim_allHLTMatch_addOfflineIso_onlyTrkIso.root'  , 'r')
# file20_iso   = TFile.Open('efficiencyReco_PF_isoMu20_TightMu_offlineRelIso_addNPrim_allHLTMatch_addOfflineIso_onlyDBIso.root'  , 'r')

c = TCanvas('', '', 600,600)
c.cd()

pt_bins  = [ 20, 24, 27, 30, 33, 36, 39, 42, 45, 49, 53, 57 ,62 ,67 ,72 ,80 , 90, 110, 150] 
eta_bins = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
 


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
    pEff1.Draw("AP")
    
    gPad.Update()

    pEff1.GetPaintedGraph().GetXaxis().SetLimits(var[5][0], var[5][1])
    pEff1.GetPaintedGraph().GetHistogram().SetMinimum(var[6][0])             
    pEff1.GetPaintedGraph().GetHistogram().SetMaximum(var[6][1])        
  
    pEff1.GetPaintedGraph().GetXaxis().SetLabelSize(0.04)
    pEff1.GetPaintedGraph().GetYaxis().SetLabelSize(0.04)   
    pEff1.GetPaintedGraph().GetXaxis().SetTitleSize(0.04)
    pEff1.GetPaintedGraph().GetYaxis().SetTitleSize(0.04)
    pEff1.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
    
    return pEff1









variables = [
#  numerator          # denominator      # x axis title                              # y axis title           # rebin      # x range      # y range  # pdf name                                    
  ('muonPt_num'       , 'muonPt_den'         , 'muon p_{T} [GeV]'                      , 'isolation efficiency', pt_bins  , (  20  ,  150), (0.8, 1),  'efficiency_muonPt'       ),
  ('muonEta_num'      , 'muonEta_den'        , 'muon #eta '                            , 'isolation efficiency', eta_bins , ( -2.4 ,  2.4), (0.8, 1),  'efficiency_muonEta'      ),
#   ('muonPhi_num'      , 'muonPhi_den'        , 'muon #phi '                            , 'isolation efficiency',       10 , ( -3.14, 3.14), (0.8, 1),  'efficiency_muonPhi_offlineRelIso0p15'      ),
#   ('muonRelIso_num'   , 'muonRelIso_den'     , 'offline #Delta#beta corrected rel iso' , 'isolation efficiency',        1 , (   0  , 0.6 ), (0. , 1),  'efficiency_muonRelIso_offlineRelIso0p15'   ),
#   ('muonRelTrkIso_num_barrel' , 'muonRelTrkIso_den_barrel'     , 'offline trk rel iso' , 'isolation efficiency',        1 , (   0  , 0.6 ), (0. , 1),  'efficiency_muonTrkRelIso_offlineRelIso0p15'   ),
#   ('muonRelPFIso_num_barrel'  , 'muonRelPFIso_den_barrel'     , 'offline ECal + HCal #Delta#beta corr. rel iso' , 'isolation efficiency',        1 , (   0  , 0.6 ), (0. , 1),  'efficiency_muonDBRelIso_offlineRelIso0p15'   ),
#   ('muonPt_num_barrel', 'muonPt_den_barrel'  , 'muon p_{T} [GeV]'                      , 'isolation efficiency', pt_bins  , (  20  ,  150), (0.8, 1),  'efficiency_muonPt_barrel_offlineRelIso0p15'),
#   ('muonPt_num_endcap', 'muonPt_den_endcap'  , 'muon p_{T} [GeV]'                      , 'isolation efficiency', pt_bins  , (  20  ,  150), (0.8, 1),  'efficiency_muonPt_endcap_offlineRelIso0p15'),
#   ('isoNvtx'          , 'muonNvtx_den'       , '# offline PV'                          , 'isolation efficiency',        1 , (  15  ,   35), (0.5, 1),  'efficiency_NPV_offlineRelIso0p15'          ),
#   ('isoNvtx_barrel'   , 'muonNvtx_den_barrel', '# offline PV'                          , 'isolation efficiency',        1 , (  15  ,   35), (0.5, 1),  'efficiency_NPV_barrel_offlineRelIso0p15'   ),
#   ('isoNvtx_endcap'   , 'muonNvtx_den_endcap', '# offline PV'                          , 'isolation efficiency',        1 , (  15  ,   35), (0.5, 1),  'efficiency_NPV_endcap_offlineRelIso0p15'   ),
]

for var in variables:
#   pEff1 = doHisto(file20, var, ROOT.kRed)
  pEff1 = doHisto(file1 , var, ROOT.kBlue)
#   pEff2 = doHisto(file20_tight, var, ROOT.kGreen+3)
#   pEff3 = doHisto(file20_iso  , var, ROOT.kOrange-3)

  c.cd()
  pEff1.Draw('AP')


  try:
    pEff2
    pEff2.Draw('P same')
    pEff3
    pEff3.Draw('P same')
    l = TLegend(0.5,0.37,0.88,0.52)
    l.SetBorderSize(0)
    l.AddEntry(pEff1 , "Reco muon + tight ID - HLT only"        , "pel")
    l.AddEntry(pEff2 , "Reco muon + tight ID - HLT only trk iso" , "pel")
    l.AddEntry(pEff3 , "Reco muon + tight ID - HLT only PFCluster iso" , "pel")
#     l.AddEntry(pEff1 , "Reco muon"        , "pel")
#     l.AddEntry(pEff2 , "Reco muon + tight ID" , "pel")
#     l.AddEntry(pEff3 , "Reco muon + tight ID + offline iso" , "pel")
    if var[0].find('barrel') > 0:  
      l.SetHeader('barrel')
    if var[0].find('endcap') > 0:  
      l.SetHeader('endcap')
    l.Draw()
  except:
    pass
    
  gPad.SetGridx(True)
  gPad.SetGridy(True)
  c.SaveAs(var[7] + "_test.pdf")







