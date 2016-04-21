# MuonHLT

cmsrel CMSSW_8_0_0  
cmsenv  
git cms-addpkg HLTrigger/Configuration  
mkdir MyTools  
cd MyTools  
git clone git@github.com:sarafiorendi/MuonHLT.git  

this gist customize the hlt configuration to run the ntuplizer 
https://gist.github.com/sarafiorendi/2e2de2da6ecbe61a218f9c203644ada2


what is in here:  
plugins/MuonNtuples.cc → code to produce ntuples   
test/evalEffArea.C     → macro to evaluate the effective areas   
test/readNtupleMC.C    → macro to evaluate efficiency on MC  
test/readNtuple.C      → macro to evaluate efficiency on data  
