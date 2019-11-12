#!/usr/bin/env python
import ROOT
import optparse

def get_histograms(files, variables):
    h = [[] for _ in range(len(files))]
    for i,f in enumerate(files):
        for var in variables:
            h[i].append( f.Get(var+'_muon;1') )
            print(h[i][-1])
    return h

def stack_histograms(h, name, legends):
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kYellow, ROOT.kGreen]
    ylow = [0.51, 0.51, 0.01, 0.01]
    yup = [0.99, 0.99, 0.49, 0.49]
    xlow = [0.01, 0.51, 0.01, 0.51]
    xup = [0.49, 0.99, 0.49, 0.99]
    assert(len(h) == len(colors)+1)

    coords1 = [0.55, 0.12, 0.12, 0.12]
    coords2 = [0.65, 0.12, 0.65, 0.65]
    coords3 = [0.71, 0.28, 0.28, 0.28]
    coords4 = [0.88, 0.35, 0.88, 0.88]
    leg = [ROOT.TLegend(coords1[i],coords2[i],coords3[i],coords4[i]) for i in range(len(h[0]))]
    c = ROOT.TCanvas('c', 'c', 1500, 1300)
    for ivar in range(len(h[0])):
        c.cd()
        pad = ROOT.TPad("pad"+str(ivar), "pad"+str(ivar), xlow[ivar], ylow[ivar], xup[ivar], yup[ivar]);
        pad.Draw()
        pad.cd()
        h[0][ivar].Draw()
        leg[ivar].AddEntry(h[0][ivar],legends[0],"f");
        for imask in range(1,len(h)):
            h[imask][ivar].SetLineColor(colors[imask-1])
            h[imask][ivar].Draw('same')
            leg[ivar].AddEntry(h[imask][ivar],legends[imask],"f");
        leg[ivar].Draw();
    c.SaveAs(name+'.png')
        

def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--outFile',          
                      dest='outFile',       
                      help='output file [%default]',  
                      default='histo_delp/val.root',       
                      type='string')
    (opt, args) = parser.parse_args()

    #inFiles = ['outMET_Mask3.root', 'outMET_Mask4.root', 'outMET_Mask5.root', 'outMET_Mask6.root', 'outMET_NoMask.root']
    inFiles = ['outMET_Mask3_v2.root', 'outMET_Mask4_v2.root', 'outMET_Mask5_v2.root', 'outMET_Mask6_v2.root', 'outMET_NoMask_v2.root']
    files = [ROOT.TFile.Open(f, 'READ') for f in inFiles]
    variables = ['MET', 'METPhi', 'Upar', 'Uperp']

    hists = get_histograms(files, variables)
    legends = ['Mask3', 'Mask4', 'Mask5', 'Mask6', 'NoMask']
    stack_histograms(hists, 'Stack', legends)

if __name__ == "__main__":
    main()
