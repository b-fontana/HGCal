#!/usr/bin/env python
import ROOT
import sys
import optparse
import numpy as np

#ROOT.gSystem.Load('../Delphes/delphes/libDelphesNoFastJet') 
ROOT.gSystem.Load('../Delphes/delphes/libDelphes') 

try:
    ROOT.gInterpreter.Declare('#include "/afs/cern.ch/user/b/bfontana/CMSSW_10_6_0/src/UserCode/Delphes/delphes/classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "/afs/cern.ch/user/b/bfontana/CMSSW_10_6_0/src/UserCode/Delphes/delphes/external/ExRootAnalysis/ExRootTreeReader.h"')
except:
    pass

detector_regions = ['End-cap',
                    'Mix',
                    'Barrel']
params = {
    "ZBosonMassRange": [76, 106],
    "plotMETRange": [-2, 25],
    "plotMETPhiRange": [-4*ROOT.TMath.Pi()/3, 4*ROOT.TMath.Pi()/3],
    "plotUparRange": [-20, 20],
    "plotUperpRange": [-20, 20],
    "plotZptRange": [0, 100],
    "plotZetaRange": [-4*ROOT.TMath.Pi()/3, 4*ROOT.TMath.Pi()/3],
    "plotEtaRange": [-5, 5],
    "plotPhiRange": [-5, 5],
    "plotMassRange": [0, 500],
    "plotNObjRange_Delp": [0, 20],
    "plotNObjRange_Full": [0, 50],
}

pdgmass = {
    'muon': 0.105658,
    'electron': 0.000511
}

nDaughters = {
    'muon': 2,
    'electron': 2
}

def createGraph(nbins, x, y, ex, ey, var):
    g = ROOT.TGraphErrors(nbins, x, y, ex, ey)
    g.GetXaxis().SetTitle( 'Region' )
    g.GetYaxis().SetTitle( 'RMS('+var+') [GeV]' )    
    g.SetMarkerColor( 1 )
    g.SetMarkerStyle( 21 )
    return g

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

def createProf(obj, varname, bins):
    prof = None
    if 'pt' in varname and 'lepton' not in varname:
        prof = ROOT.TProfile(varname+'_'+obj, "", len(bins)-1, bins)
        prof.GetXaxis().SetTitle("p_{T}(Z) [GeV]")
        prof.GetYaxis().SetTitle('<'+varname[3:]+'>')
    elif 'eta' in varname and 'lepton' not in varname:
        prof = ROOT.TProfile(varname+'_'+obj, "", len(bins)-1, bins)
        prof.GetXaxis().SetTitle("#eta(Z)")
        prof.GetYaxis().SetTitle('<'+varname[4:]+'>')
    elif 'lepton_eta' in varname:
        prof = ROOT.TProfile(varname+'_'+obj, "", len(bins)-1, bins)
        prof.GetXaxis().SetTitle("Region")
        prof.GetXaxis().SetBinLabel(1, detector_regions[0])
        prof.GetXaxis().SetBinLabel(2, detector_regions[1])
        prof.GetXaxis().SetBinLabel(3, detector_regions[2])
        prof.GetXaxis().CenterLabels()
        prof.GetYaxis().SetTitle('<'+varname[11:]+'>')
    if prof==None:
        raise ValueError('No profile was created.')
    return prof

"""
def check_size(event):
    sizes = {'muon': 2}
    for k,v in sizes.items():
        if k == 'muon':
            if event.Muon_size != v:
                return 1
        else:
            raise ValueError('Muons only!')
    return 0
"""

def main():
    objects = ['muon', 'electron'] 
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
                      default=999999,
                      type=int)
    (opt, args) = parser.parse_args()

    inFile = opt.inFile    
    f = ROOT.TFile.Open(inFile, 'READ')
    ntuple = f.Get('Delphes;1')

    maxEvents = opt.maxEvts

    tot_nevents = 0
    tot_muons = 0
    tot_electrons = 0

    muon_mothers = []

    outputF = ROOT.TFile(opt.outFile, "RECREATE")
    hists = {} 
    profs = {}
    rms_hists = {}
    pt_bins = np.array([0.,5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,100.])
    eta_bins = np.array([-8*ROOT.TMath.Pi()/6, -7*ROOT.TMath.Pi()/6, -ROOT.TMath.Pi(), -5*ROOT.TMath.Pi()/6, -4*ROOT.TMath.Pi()/6, -ROOT.TMath.Pi()/2, -2*ROOT.TMath.Pi()/6, -ROOT.TMath.Pi()/6, 0., ROOT.TMath.Pi()/6, 2*ROOT.TMath.Pi()/6, ROOT.TMath.Pi()/2, 4*ROOT.TMath.Pi()/6, 5*ROOT.TMath.Pi()/6, ROOT.TMath.Pi(), 7*ROOT.TMath.Pi()/6, 8*ROOT.TMath.Pi()/6])
    lepton_eta_bins = np.array([0.,1.,2.,3.])
    for hname in variables:
        for obj in objects:
            hists[hname+'_'+obj] = createHist(obj, hname)
            profs['pt'+'_'+hname+'_'+obj] = createProf(obj, 'pt'+'_'+hname, pt_bins)
            profs['eta'+'_'+hname+'_'+obj] = createProf(obj, 'eta'+'_'+hname, eta_bins)
            profs['lepton_eta'+'_'+hname+'_'+obj] = createProf(obj, 'lepton_eta'+'_'+hname, lepton_eta_bins)
            if hname == 'Upar' or hname == 'Uperp':
                for irms in range(3):
                    rms_hists['rms_'+hname+'_'+obj+str(irms)] = createHist(obj, hname)

    for entry,event in enumerate(ntuple):
        if maxEvents > 0 and entry >= maxEvents:
            break
        if (tot_nevents %100) == 0 :
            print('... processed {} events ...'.format(entry+1))

        tot_nevents += 1
        
        """
        if check_size(event) == 1:
            continue
        """

        ######################################
        lepton_barrel = dict()
        for obj in objects:
            lepton_barrel.update({obj: 0}) #initialize to zero
        mother = dict()

        for obj in objects:
            if obj == 'muon':
                event_obj = event.Muon
                tot_muons += event.Muon_size
            elif obj == 'electron':
                event_obj = event.Electron
                tot_electrons += event.Electron_size
            else:
                raise ValueError('Muons and eletrons only!')

            obj_vectors = []        
            for p in event_obj:
                if (abs(p.Eta) < 1.45):
                    lepton_barrel[obj] += 1
                obj_vectors.append( ROOT.Math.PtEtaPhiMVector(p.PT, p.Eta, p.Phi, pdgmass[obj]) )

            #reconstruct mother particle if the number of daughters is correct
            if len(obj_vectors) == nDaughters[obj]:
                mother[obj] = obj_vectors[0]
                for vector in obj_vectors[1:]:
                    mother[obj] += vector

        ######################################

        for obj in objects:
            if obj not in mother: #the mother was not reconstructed for this object
                continue
            if mother[obj].M() > params["ZBosonMassRange"][0] and mother[obj].M() < params["ZBosonMassRange"][1]:
                for ev in event.MissingET: #one single value per event
                    phi_par = ROOT.TVector2.Phi_mpi_pi( ev.Phi - (mother[obj].Phi() + ROOT.TMath.Pi()) )
                    upar = ev.MET * ROOT.TMath.Cos(phi_par)
                    uperp = ev.MET * ROOT.TMath.Sin(phi_par)

                    hists['MET_'+obj].Fill(ev.MET)
                    hists['METPhi_'+obj].Fill(ev.Phi)
                    hists['Upar_'+obj].Fill(upar)
                    hists['Uperp_'+obj].Fill(uperp)
                    
                    profs['pt_MET_'+obj].Fill(mother[obj].Pt(), ev.MET)
                    profs['pt_METPhi_'+obj].Fill(mother[obj].Pt(), ev.Phi)
                    profs['pt_Upar_'+obj].Fill(mother[obj].Pt(), upar)
                    profs['pt_Uperp_'+obj].Fill(mother[obj].Pt(), uperp)

                    profs['eta_MET_'+obj].Fill(mother[obj].Eta(), ev.MET)
                    profs['eta_METPhi_'+obj].Fill(mother[obj].Eta(), ev.Phi)
                    profs['eta_Upar_'+obj].Fill(mother[obj].Eta(), upar)
                    profs['eta_Uperp_'+obj].Fill(mother[obj].Eta(), uperp)

                    profs['lepton_eta_MET_'+obj].Fill(lepton_barrel[obj], ev.MET)
                    profs['lepton_eta_METPhi_'+obj].Fill(lepton_barrel[obj], ev.Phi)
                    profs['lepton_eta_Upar_'+obj].Fill(lepton_barrel[obj], upar)
                    profs['lepton_eta_Uperp_'+obj].Fill(lepton_barrel[obj], uperp)

                    rms_hists['rms_Upar_'+obj+str(lepton_barrel[obj])].Fill(upar) 
                    rms_hists['rms_Uperp_'+obj+str(lepton_barrel[obj])].Fill(uperp) 

    
    
    #Extract the RMS. Each histogram corresponds to one 'detector_regions' bin
    for obj in objects:
        for var in ['Upar', 'Uperp']:
            arr = np.array([rms_hists['rms_'+var+'_'+obj+'0'].GetRMS(),
                            rms_hists['rms_'+var+'_'+obj+'1'].GetRMS(),
                            rms_hists['rms_'+var+'_'+obj+'2'].GetRMS()])
            arr_e = np.array([rms_hists['rms_'+var+'_'+obj+'0'].GetRMSError(),
                              rms_hists['rms_'+var+'_'+obj+'1'].GetRMSError(),
                              rms_hists['rms_'+var+'_'+obj+'2'].GetRMSError()])
            g = createGraph(3, np.array([1,2,3], dtype='d'), arr, 
                            np.array([0.5, 0.5, 0.5]), arr_e, var)
            g.Draw( 'AP' )
            g.SetName('rms_'+var+'_'+obj)
            outputF.cd()
            g.Write()

    #Save
    outputF.cd()
    for h in hists.keys():
        hists[h].Write()
    for p in profs.keys():
        profs[p].Write()

    print("Processed %d events" % tot_nevents)
    print("On average %f muons" % (float(tot_muons) / tot_nevents))
    print("On average %f electrons" % (float(tot_electrons) / tot_nevents))

if __name__ == "__main__":
    main()
