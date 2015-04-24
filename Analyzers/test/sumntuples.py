import os
import sys 
import ROOT
import copy
import glob
from ROOT import TFile, gDirectory
import subprocess

def listfiles(eosfolder) :
  proc = subprocess.Popen(['/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select','ls',eosfolder],stdout=subprocess.PIPE)
  files = []
  for line in proc.stdout:
    files.append('root://eoscms//eos/cms/'+eosfolder+'/'+line.rstrip())
  return files

if len(sys.argv) != 3 :
  print 'usage:     python sumMCFiles.py folder outfile'
  exit(1)

toSum     = []
eosfolder = sys.argv[1]
outfile   = sys.argv[2]
files = listfiles(eosfolder)
for file in files : 
#   print file
  if 'ntuple' in file :
    toSum.append(file)

#for ifile in toSum:
#  print ifile

#exit(1)
cmd1 = 'hadd ' + outfile + ' ' 
for ifile in toSum:
  cmd1 += str(ifile)+ ' '

print cmd1
os.system(cmd1)
print ' '
# 
# 
# if len(argv) != 3 :
#  print 'usage:     python sumMCFiles.py folder outfile'
#  exit(1)
# 
# 
# folder  = argv[1]
# outfile = argv[2]
# print outfile
# 
# 
# toSum     = []
# counter = 0
# dirList = os.listdir(folder)
# for fname in dirList:
#   print fname
#   if  os.path.isdir(fname) : continue
#   if 'ntuple' in fname :
#     toSum.append(folder+'/'+fname)


# print 'I\'m going to add histos for folder: ', str(folder)
# 	  
# cmd1 = 'hadd ' + folder + '/treeMC' + sample + '.root '
# for ifile in toSum:
#   cmd1 += str(ifile)+ ' '
# # print cmd1
# 
# # os.system(cmd1)
# print ' '
# print ''
