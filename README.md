# MuonHLT

```
cmsrel CMSSW_10_1_7  
cd CMSSW_10_1_7/src  
cmsenv    
git cms-addpkg HLTrigger/Configuration    
git clone git@github.com:sarafiorendi/MuonHLT.git    
cd MuonHLT/  
git checkout KeeStudies   
cd ..  
scramv1 b   
```

To produce ntuples, just run    
`cmsRun hltNtuples_cfg.py`
from the test folder

To re-run the hlt and produce the ntuples:
```
cd $CMSSW_BASE/src/HLTrigger/Configuration/test
hltGetConfiguration /dev/CMSSW_10_1_0/GRun .... > hlt_parking.py
cmsRun customize_hlt_purity.py
```

where customize_hlt_purity.py is something similar to this gist:
https://gist.github.com/sarafiorendi/6b6b6f6702265eb94bb2b0e4f80c014d
