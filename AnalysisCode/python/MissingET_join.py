#!/usr/bin/env python
import os
import ROOT
from ROOT import Double
import optparse
import copy
stats_file = 'stats.txt'

def get_histograms(files, variables, obj):
    h = [[] for _ in range(len(files))]
    for i,f in enumerate(files):
        for var in variables:
            h[i].append( f.Get(var+'_'+obj+';1') )
            print(h[i][-1])
    return h

def print_stats(h, mask_str):
    assert(len(h) == 3)
    mean = [round(k.GetMean(), 3) for k in h]
    emean = [round(k.GetMeanError(), 3) for k in h]
    rms = [round(k.GetRMS(), 3) for k in h]
    erms = [round(k.GetRMSError(), 3) for k in h]
    with open(stats_file, 'a') as f:                                                                                                                                         
        f.write('{} & $\mu = {} \pm {}$ & $\mu = {} \pm {}$ & $\mu = {} \pm {}$ \\\ \n'
                .format(mask_str,mean[0],emean[0],mean[1],emean[1],mean[2],emean[2]))
        f.write(' & $\sigma = {} \pm {}$ & $\sigma = {} \pm {}$ & $\sigma = {} \pm {}$ \\\ \n'
                .format(rms[0],erms[0],rms[1],erms[1],rms[2],erms[2]))

def stack_histograms(h, name, legends):
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kYellow, ROOT.kGreen]
    assert(len(h) == len(colors)+1)
    minimums, maximums = [-1, -3.9, -19, -19], [24.5, 3.9, 19.5, 19.5]
    coords = [[0.55, 0.12, 0.12, 0.12],
              [0.65, 0.12, 0.65, 0.65],
              [0.71, 0.28, 0.28, 0.28],
              [0.88, 0.35, 0.88, 0.88]]
    leg = [ROOT.TLegend(coords[0][i],coords[1][i],coords[2][i],coords[3][i]) for i in range(len(h[0]))]

    for ivar in range(len(h[0])):
        for imask in range(0,len(h)):
            #h[imask][ivar].SetAxisRange(minimums[ivar]-10, maximums[ivar]+10,"X");
            h[imask][ivar].GetXaxis().SetRangeUser(0, 2);
            h[imask][ivar] = h[imask][ivar].Rebin(2)
            print(ivar, imask, h[imask][ivar].GetEntries())
            scale = 1. / h[imask][ivar].Integral()
            h[imask][ivar].Scale(scale)

    for imask in range(0,len(h)):
        print_stats([h[imask][0],h[imask][2],h[imask][3]], legends[imask])
            
    for ivar in range(len(h[0])):
        c = ROOT.TCanvas('c'+str(ivar), 'c'+str(ivar), 1500, 1300)
        c.cd()
        pad = ROOT.TPad('pad'+str(ivar), 'pad'+str(ivar), 0., 0.36, 1., 1.);
        pad.SetBottomMargin(0.);
        pad.Draw()
        pad.cd()
        h[0][ivar].Draw()
        leg[ivar].AddEntry(h[0][ivar],legends[0],"f");
        for imask in range(1,len(h)):
            h[imask][ivar].SetLineColor(colors[imask-1])
            h[imask][ivar].Draw('same')
            h[imask][ivar].GetXaxis().SetLabelSize(0);
            leg[ivar].AddEntry(h[imask][ivar],legends[imask],"f");
        leg[ivar].Draw();
        c.cd()
        pad2 = ROOT.TPad('pad'+str(ivar)+'_low', 'pad'+str(ivar)+'_low', 0., 0.01, 1., 0.36);
        pad2.Draw()
        pad2.cd()
        h_clone = []
        for imask in range(1,len(h)):
            h_clone.append(h[imask][ivar].Clone('h_clone_mask'+str(imask+2)+'_var'+str(ivar)))
            h_clone[-1].Divide(h[0][ivar])
            h_clone[-1].SetLineColor(colors[imask-1])
            h_clone[-1].SetStats(ROOT.kFALSE)
            ax, ay = h_clone[-1].GetXaxis(), h_clone[-1].GetYaxis()
            ay.SetRangeUser(0.5,1.5);
            ay.SetTitle("");
            ax.SetLabelSize(0.05);
            ax.SetTitleSize(0.05);
            ax.SetTitleOffset(0.95);
            ay.SetLabelSize(0.06);
            ay.SetNdivisions(-202);
            h_options = 'l'
            h_clone[-1].Draw('hist same '+h_options) if imask!=1 else h_clone[-1].Draw('hist '+h_options)
        l = ROOT.TLine()
        l.SetLineStyle(9)
        l.SetLineWidth(4)
        l.SetLineColor(15) #grey
        l.DrawLine(minimums[ivar], 1., maximums[ivar], 1.)
        c.cd()
        pad.cd()
        redraw_border()
        c.SaveAs('/eos/user/b/bfontana/www/PartialWafers/MissingET/'+name+str(ivar)+'.png')
        c.Close() #needed to avoid repeated name error

def stack_profiles(g, name, legends, minimum, maximum):
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kYellow, ROOT.kGreen]
    assert(len(g) == len(colors)+1)
    if 'pt' in name and 'lepton' not in name:
        coords = [[0.73, 0.12, 0.12, 0.12],
                  [0.01, 0.12, 0.65, 0.65],
                  [0.89, 0.28, 0.28, 0.28],
                  [0.24, 0.35, 0.88, 0.88]]
    elif 'eta' in name and 'lepton' not in name:
        coords = [[0.73, 0.12, 0.12, 0.12],
                  [0.65, 0.12, 0.65, 0.65],
                  [0.89, 0.28, 0.28, 0.28],
                  [0.88, 0.35, 0.88, 0.88]]
    else:
        coords = [[0.73, 0.12, 0.73, 0.73],
                  [0.76, 0.12, 0.76, 0.76],
                  [0.89, 0.28, 0.89, 0.89],
                  [0.89, 0.35, 0.89, 0.89]]
    leg = [ROOT.TLegend(coords[0][i],coords[1][i],coords[2][i],coords[3][i]) for i in range(len(g[0]))]

    for ivar in range(len(g[0])):
        c = ROOT.TCanvas('c'+str(ivar), 'c'+str(ivar), 1500, 1300)
        c.cd()
        pad = ROOT.TPad('pad'+str(ivar), 'pad'+str(ivar), 0., 0.36, 1., 1.);
        pad.SetBottomMargin(0.);
        pad.Draw()
        pad.cd()
        if 'lepton_eta' in name:
            define_lepton_limits(g[0][ivar], ivar)
        g[0][ivar].Draw()
        leg[ivar].AddEntry(g[0][ivar],legends[0],"f");
        for imask in range(1,len(g)):
            g[imask][ivar].SetLineColor(colors[imask-1])
            g[imask][ivar].Draw('same')
            g[imask][ivar].GetXaxis().SetLabelSize(0);
            leg[ivar].AddEntry(g[imask][ivar],legends[imask],"f");
        leg[ivar].Draw();
        c.cd()
        pad2 = ROOT.TPad('pad'+str(ivar)+'_low', 'pad'+str(ivar)+'_low', 0., 0.01, 1., 0.36);
        pad2.Draw()
        pad2.cd()
        g_clone = []
        for imask in range(1,len(g)):
            g_clone.append(g[imask][ivar].Clone('g_clone_mask'+str(imask+2)+'_var'+str(ivar)))
            g_clone[-1].Divide(g[0][ivar])
            g_clone[-1].SetLineColor(colors[imask-1])
            g_clone[-1].SetStats(ROOT.kFALSE)
            ax, ay = g_clone[-1].GetXaxis(), g_clone[-1].GetYaxis()
            ay.SetRangeUser(0.94,1.06);
            ay.SetTitle("");
            ax.SetLabelSize(0.05);
            ax.SetTitleSize(0.05);
            ax.SetTitleOffset(0.95);
            ay.SetLabelSize(0.06);
            ay.SetNdivisions(-202);
            g_options = 'l'
            g_clone[-1].Draw('hist same '+g_options) if imask!=1 else g_clone[-1].Draw('hist '+g_options)
        l = ROOT.TLine()
        l.SetLineStyle(9)
        l.SetLineWidth(4)
        l.SetLineColor(15) #grey
        l.DrawLine(minimum, 1., maximum, 1.)
        c.cd()
        pad.cd()
        redraw_border()
        c.SaveAs('/eos/user/b/bfontana/www/PartialWafers/MissingET/'+name+str(ivar)+'.png')
        c.Close() #needed to avoid repeated name error

def stack_graphs(g, name, legends, minimum, maximum):
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kYellow, ROOT.kGreen]
    assert(len(g) == len(colors)+1)
    coords = [[0.73, 0.73],
              [0.76, 0.76],
              [0.89, 0.89],
              [0.89, 0.89]]
    leg = [ROOT.TLegend(coords[0][i],coords[1][i],coords[2][i],coords[3][i]) for i in range(len(g[0]))]
     
    for ivar in range(len(g[0])):
        c = ROOT.TCanvas('c'+str(ivar), 'c'+str(ivar), 1500, 1300)
        c.cd()
        pad = ROOT.TPad('pad'+str(ivar), 'pad'+str(ivar), 0., 0.36, 1., 1.);
        pad.SetBottomMargin(0.);
        pad.Draw()
        pad.cd()
        define_rms_limits(g[0][ivar], ivar)
        g[0][ivar].Draw()
        leg[ivar].AddEntry(g[0][ivar],legends[0],"f");
        for imask in range(1,len(g)):
            g[imask][ivar].SetLineColor(colors[imask-1])
            g[imask][ivar].Draw('same')
            g[imask][ivar].GetXaxis().SetLabelSize(0);
            leg[ivar].AddEntry(g[imask][ivar],legends[imask],"f");
        leg[ivar].Draw();
        c.cd()
        pad2 = ROOT.TPad('pad'+str(ivar)+'_low', 'pad'+str(ivar)+'_low', 0., 0.01, 1., 0.36);
        pad2.Draw()
        pad2.cd()
        g_clone = [ROOT.TGraphErrors(g[0][ivar].GetN()) for _ in range(len(g))]
        for imask in range(1,len(g)):
            for ipoint in range(g[0][ivar].GetN()):
                pointx, pointx0 = Double(0.), Double(0.)
                pointy, pointy0 = Double(0.), Double(0.)
                errorx, errorx0 = Double(0.), Double(0.)
                errory, errory0 = Double(0.), Double(0.)
                g[0][ivar].GetPoint(ipoint, pointx0, pointy0)
                g[imask][ivar].GetPoint(ipoint, pointx, pointy)
                assert(pointx == pointx0)
                g_clone[imask-1].SetPoint(ipoint, pointx, pointy/pointy0)
                errorx0 = g[0][ivar].GetErrorX(ipoint)
                errory0 = g[0][ivar].GetErrorY(ipoint)
                errorx = g[imask][ivar].GetErrorX(ipoint)
                errory = g[imask][ivar].GetErrorY(ipoint)
                assert(errorx == errorx0)
                newerrory = (pointy/pointy0) * ROOT.TMath.Sqrt( (errory/pointy)**2 + (errory0/pointy0)**2 )
                g_clone[imask-1].SetPointError(ipoint, errorx, newerrory)
            g_clone[imask-1].SetMarkerSize(1.5)
            g_clone[imask-1].SetMarkerColor(colors[imask-1])
            g_clone[imask-1].SetLineColor(colors[imask-1])
            g_clone[imask-1].SetMarkerStyle(8)
            ax, ay = g_clone[imask-1].GetXaxis(), g_clone[imask-1].GetYaxis()
            ay.SetRangeUser(0.94,1.06);
            ay.SetTitle("");
            ax.SetLabelSize(0.05);
            ax.SetTitleSize(0.05);
            ax.SetTitleOffset(0.95);
            ay.SetLabelSize(0.06);
            ay.SetNdivisions(-202);
            g_options = 'p'
            g_clone[imask-1].Draw('same '+g_options) if imask!=1 else g_clone[imask-1].Draw('a'+g_options)
        l = ROOT.TLine()
        l.SetLineStyle(9)
        l.SetLineWidth(4)
        l.SetLineColor(15) #grey
        l.DrawLine(minimum, 1., maximum, 1.)
        c.cd()
        pad.cd()
        redraw_border()
        c.SaveAs('/eos/user/b/bfontana/www/PartialWafers/MissingET/'+name+str(ivar)+'.png')
        c.Close() #needed to avoid repeated name error
        
def redraw_border():
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    l = ROOT.TLine()
    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    l2 = ROOT.TLine()
    l2.DrawLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())

def define_lepton_limits(graph, i):
    if i==0:
        graph.GetYaxis().SetRangeUser(5.6, 8.)
    elif i==1:
        graph.GetYaxis().SetRangeUser(5.6, 8.)
    elif i==2:
        graph.GetYaxis().SetRangeUser(-3.3, 3.3)
    elif i==3:
        graph.GetYaxis().SetRangeUser(-0.3, 0.3)
    else:
        raise ValueError('Wrong variable number.')

def define_rms_limits(graph, i):
    if i==0:
        graph.GetYaxis().SetRangeUser(4.8, 6.35)
    elif i==1:
        graph.GetYaxis().SetRangeUser(4.4, 5.5)
    else:
        raise ValueError('Wrong variable number.')

def main():
    ROOT.gStyle.SetOptStat(0);
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    try:
        os.remove(stats_file)
    except OSError:
        pass

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--outFile',          
                      dest='outFile',       
                      help='output file [%default]',  
                      default='histo_delp/val.root',       
                      type='string')
    (opt, args) = parser.parse_args()

    objects = ['muon', 'electron']
    legends = ['NoMask', 'Mask3', 'Mask4', 'Mask5', 'Mask6']
    inFiles = ['outMET_'+leg+'_ee_mumu.root' for leg in legends]
    files = [ROOT.TFile.Open(f, 'READ') for f in inFiles]
    variables = ['MET', 'METPhi', 'Upar', 'Uperp']
    variables_pt = ['pt_MET', 'pt_METPhi', 'pt_Upar', 'pt_Uperp']
    variables_eta = ['eta_MET', 'eta_METPhi', 'eta_Upar', 'eta_Uperp']
    variables_leptons = ['lepton_eta_MET', 'lepton_eta_METPhi', 'lepton_eta_Upar', 'lepton_eta_Uperp']
    variables_rms = ['rms_Upar', 'rms_Uperp']

    for obj in objects:
        hists = get_histograms(files, variables, obj)
        stack_histograms(hists, obj+'_', legends)
        profiles_pt = get_histograms(files, variables_pt, obj)
        stack_profiles(profiles_pt, 'pt_'+obj+'_', legends, -2, 122)
        profiles_eta = get_histograms(files, variables_eta, obj)
        stack_profiles(profiles_eta, 'eta_'+obj+'_', legends, -4.5, 4.5)
        profiles_leptons = get_histograms(files, variables_leptons, obj)
        stack_profiles(profiles_leptons, 'lepton_eta_'+obj+'_', legends, 0.2, 2.8)
        graphs_rms = get_histograms(files, variables_rms, obj)
        stack_graphs(graphs_rms, 'rms_'+obj+'_', legends, 0.5, 3.6)

if __name__ == "__main__":
    main()
