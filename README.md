# MuonHLT

cmsrel CMSSW_8_0_8_patch1  
cd CMSSW_8_0_8_patch1/src  
cmsenv    
git cms-addpkg HLTrigger/Configuration    
git clone git@github.com:sarafiorendi/MuonHLT.git    
cd MuonHLT/  
git checkout iterL3 
#git merge origin/addL1info  
cd ..  
scramv1 b   
  
to produce ntuples, just run    
cmsRun hltNtuples_cfg.py

this gist customize the hlt configuration to re-run the hlt and the ntuplizer for isolation studies
https://gist.github.com/sarafiorendi/02b4a43c28766ec7beddabe3f5ae036e
  
what is in here:  
plugins/MuonNtuples.cc → code to produce ntuples   
test/evalEffArea.C     → macro to evaluate the effective areas on MC   
test/readNtupleMC.C    → macro to evaluate efficiency on MC  
test/readNtuple.C      → macro to evaluate efficiency on data  
