import ROOT
import sys
import os
inDir=sys.argv[1]
allFiles=[os.path.join(inDir,f) for f in os.listdir(inDir) if '.root' in f]
for f in allFiles:
    inF=ROOT.TFile.Open(f)
    toDel=False
    try:
        if inF.IsZombie(): toDel=True
        nevts=inF.Get('Events').GetEntriesFast()
        if nevts==0: toDel=True
    except:
        toDel=True
    if not toDel: continue    
    os.system('rm -v %s'%f)
