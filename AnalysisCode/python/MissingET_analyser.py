#!/usr/bin/env python
import ROOT
import sys
import optparse

#ROOT.gSystem.Load('../Delphes/delphes/libDelphesNoFastJet') 
ROOT.gSystem.Load('../Delphes/delphes/libDelphes') 

try:
    ROOT.gInterpreter.Declare('#include "/afs/cern.ch/user/b/bfontana/CMSSW_10_6_0/src/UserCode/Delphes/delphes/classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "/afs/cern.ch/user/b/bfontana/CMSSW_10_6_0/src/UserCode/Delphes/delphes/external/ExRootAnalysis/ExRootTreeReader.h"')
except:
    pass

params = {
    "ZBosonMassRange": [106, 76],
    "plotMETRange": [-2, 25],
    "plotMETPhiRange": [-4*ROOT.TMath.Pi()/3, 4*ROOT.TMath.Pi()/3],
    "plotUparRange": [-20, 20],
    "plotUperpRange": [-20, 20],
    "plotEtaRange": [-5, 5],
    "plotPhiRange": [-5, 5],
    "plotMassRange": [0, 500],
    "plotNObjRange_Delp": [0, 20],
    "plotNObjRange_Full": [0, 50],
}

pdgmass = {
    'muon': 0.105658,
}

def createHist(obj, varname):
    h = None
    if varname == 'MET':
        h = ROOT.TH1D(varname+'_'+obj, "", 100, params["plotMETRange"][0], params["plotMETRange"][1])
        h.GetXaxis().SetTitle("MET [GeV]")
        h.GetYaxis().SetTitle("N")
    elif varname == 'METPhi':
        h = ROOT.TH1D(varname+'_'+obj, "", 100, params["plotMETPhiRange"][0], params["plotMETPhiRange"][1])
        h.GetXaxis().SetTitle("METPhi [GeV]")
        h.GetYaxis().SetTitle("N")
    elif varname == 'Upar':
        h = ROOT.TH1D(varname+'_'+obj, "", 100, params["plotUparRange"][0], params["plotUparRange"][1])
        h.GetXaxis().SetTitle("Upar [GeV]")
        h.GetYaxis().SetTitle("N")
    elif varname == 'Uperp':
        h = ROOT.TH1D(varname+'_'+obj, "", 100, params["plotUperpRange"][0], params["plotUperpRange"][1])
        h.GetXaxis().SetTitle("Uperp [GeV]")
        h.GetYaxis().SetTitle("N")
    if h==None:
        raise ValueError('No histogram was created.')
    h.Sumw2()
    return h

def check_size(event):
    sizes = {'muon': 2}
    for k,v in sizes.items():
        if k == 'muon':
            if event.Muon_size != v:
                return 1
        else:
            raise ValueError('Muons only!')
    return 0

def main():
    objects = ['muon'] 
    variables = ['MET', 'METPhi', 'Upar', 'Uperp']

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inFile',
                      dest='inFile',
                      help='input file [%default]',  
                      default=None,
                      type='string')
    parser.add_option('-o', '--outFile',          
                      dest='outFile',       
                      help='output file [%default]',  
                      default='histo_delp/val.root',       
                      type='string')
    parser.add_option('--maxEvents',          
                      dest='maxEvts',
                      help='max number of events [%default]',
                      default=1000,
                      type=int)
    (opt, args) = parser.parse_args()

    inFile = opt.inFile    
    f = ROOT.TFile.Open(inFile, 'READ')
    ntuple = f.Get('Delphes;1')

    maxEvents = opt.maxEvts

    tot_nevents = 0
    tot_muons = 0

    muon_mothers = []

    outputF = ROOT.TFile(opt.outFile, "RECREATE")
    hists = {} 
    for hname in variables:
        for obj in objects:
            hists[hname] = createHist(obj, hname)

    for entry,event in enumerate(ntuple):
        if maxEvents > 0 and entry >= maxEvents:
            break
        if (tot_nevents %100) == 0 :
            print('... processed {} events ...'.format(entry+1))

        tot_nevents += 1

        if check_size(event) == 1:
            continue

        ######################################
        #for/else for applying cuts#
        mother = dict()
        for obj in objects:
            if obj == 'muon':
                event_obj = event.Muon
                tot_muons += event.Muon_size
            else:
                raise ValueError('Muons only!')

            obj_vectors = []        
            for p in event_obj:
                obj_vectors.append( ROOT.Math.PtEtaPhiMVector(p.PT, p.Eta, p.Phi, pdgmass[obj]) )
            if len(obj_vectors) != 2:
                break;

            mother[obj] = obj_vectors[0]
            for vector in obj_vectors[1:]:
                mother[obj] += vector

            if mother[obj].M() < params["ZBosonMassRange"][0] or mother[obj].M() > params["ZBosonMassRange"][1]:
                break

        else:
            continue
        ######################################

        #one single value per event
        for ev in event.MissingET:
            phi_par = ROOT.TVector2.Phi_mpi_pi( ev.Phi - (mother['muon'].Phi() + ROOT.TMath.Pi()) )
            hists['MET'].Fill(ev.MET)
            hists['METPhi'].Fill(ev.Phi)
            hists['Upar'].Fill(ev.MET * ROOT.TMath.Cos(phi_par))
            hists['Uperp'].Fill(ev.MET * ROOT.TMath.Sin(phi_par))
            
        # end of one event

    outputF.cd()
    for h in hists.keys():
        hists[h].Write()

    print("Processed %d events" % tot_nevents)
    print("On average %f muons" % (float(tot_muons) / tot_nevents))

if __name__ == "__main__":
    main()
