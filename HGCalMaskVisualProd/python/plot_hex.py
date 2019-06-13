import sys, os
import numpy as np
import argparse
import math, copy
from matplotlib import pyplot as plt

from ROOT import TCanvas, TFile, TIter, TH2Poly, TPaveLabel, TKey
from ROOT import gDirectory, gROOT, gStyle, kTRUE, kFALSE
gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(kFALSE)
gStyle.SetPalette(53)
gStyle.SetNumberContours(99)

import RootUtils, SystemUtils, Argparser

parser = argparse.ArgumentParser()
FLAGS, _ = Argparser.add_args(parser)
Argparser.print_args(FLAGS)

def getUVFromTitle(title):
    u, v, N, rest = title.split(',')
    return u, v, N, rest

def convertToHex_alternative(h, N, a):
    """
    Arguments:
    -> h: the standard histogram
    -> N: N as defined by the HGCal geometry
    -> a: the length in cm of each side of each cell
    """
    def convertCoords(u, v):
        return 1.5*(v-N)+1, u-0.5*(N+v)

    def allowedWafers():
        m = []
        for u,v in zip(range(N+1),range(N-1,2*N,1)):
            for i in range(v+1):
                m.append(tuple((u,i)))
        for u,v in zip(range(N+1,2*N),range(1,N)):
            for i in range(v,2*N):
                m.append(tuple((u,i)))
        return m
        

    def addBin(xc, yc):
        """
        Arguments:
        -> xc: x coordinate of the centre of the bin
        -> yc: y coordinate of the centre of the bin
        """
        from array import array
        x = array('d', [0., 0., 0., 0., 0., 0.])
        y = array('d', [0., 0., 0., 0., 0., 0.])
        x[0] = xc - a
        y[0] = yc
        x[1] = x[0] + a/2.
        y[1] = y[0] + a*math.sqrt(3)/2.
        x[2] = x[1] + a
        y[2] = y[1]
        x[3] = x[0] + 2*a
        y[3] = y[0]
        x[4] = x[2]
        y[4] = y[3] - a*math.sqrt(3)/2.
        x[5] = x[1]
        y[5] = y[4] 
        hhex.AddBin(6, x, y)
 
    
    lim = 20 if N==12 else 13
    hhex = TH2Poly('hhex', 'hhex', -lim, lim, -lim, lim)
    xylist = []
    allWaf = allowedWafers()

    for U in range(h.GetNbinsX()):
        for V in range(h.GetNbinsY()):
            if (U,V) in allWaf:
                x, y = convertCoords(U, V)
                if (x,y) not in xylist:
                    addBin(x, y)
                    xylist.append((x,y))

    for U in range(h.GetNbinsX()):
        for V in range(h.GetNbinsY()):
            x, y = convertCoords(U, V)
            if (x,y) in xylist:  
                global_bin_standard = h.FindBin(U, V)
                content_standard = h.GetBinContent(global_bin_standard)
                global_bin_hex = hhex.FindBin(x, y)
                hhex.SetBinContent(global_bin_hex, content_standard)
    return hhex

def keyTitleMap(keylist):
    """
    Creates a map linking histogram titles to keys.
    Arguments:
    -> keylist: The list of keys from which the map is created
    Returns:
    -> a map
    """
    m = {}
    for k in keylist:
        if not gROOT.GetClass(k.GetClassName()).InheritsFrom('TH1'):
            continue
        m.update({k.ReadObj().GetTitle(): k})
    return m


pcoords = [[[0.02,0.49,0.98,0.93],    #pad1
           [0.02,0.02,0.98,0.48]]]   #pad2
pad_title = ['RecHits', 'Geom']
for ic,_ in enumerate(pcoords):
    if len(pcoords[ic]) != len(pad_title):
        print(ic, len(pcoords[ic]), len(pad_title))
        raise ValueError('The number of pads has to be consistent.')

myfile = TFile.Open(FLAGS.root_file, "READ")
online_dir = '/eos/user/b/bfontana/www/'
maind = 'outpath'
layers = [x for x in range(1,FLAGS.nlayers+1)]
subd_names = ['layer'+str(layers[x]) for x in range(len(layers))]
for isd,sd in enumerate(subd_names):
    directory = os.path.join(online_dir,FLAGS.directory,sd)
    SystemUtils.createDir(os.path.join(online_dir,FLAGS.directory))
    SystemUtils.createDir(directory)
    myfile.cd(os.path.join(maind,sd))
    kmap = keyTitleMap(gDirectory.GetListOfKeys())
    titleUsed = []
    for htitle,key in kmap.items():
        h = key.ReadObj()
        uvN = getUVFromTitle(htitle)
        title1 = uvN[0]+','+uvN[1]+','+uvN[2]+',RecHits'
        title2 = uvN[0]+','+uvN[1]+','+uvN[2]+',Geom'
        if title1 in titleUsed or title2 in titleUsed:
            continue
        titleUsed.append([title1, title2])
        if htitle==title1:
            histos = [h, kmap[title2].ReadObj()]
        elif htitle==title2:
            histos = [kmap[title1].ReadObj(), h]
        else:
            raise ValueError('There is a problem with the title of the histogram.')

        with RootUtils.RootPlotting(ncanvas=1, npads=len(pad_title), 
                                    pcoords=pcoords) as plot:
            #for ipad in range(len(pad_title)):
            ipad=0
            h1 = convertToHex_alternative(histos[0], int(uvN[2]), 0.5)
            plot.plotHistogram(cpos=0, ppos=ipad, h=h1, title=pad_title[ipad], 
                               xaxis_title='X [cm]', yaxis_title='Y [cm]',
                               draw_options='colz')

            ipad=1
            h2 = convertToHex_alternative(histos[1], int(uvN[2]), 0.5)
            plot.plotHistogram(cpos=0, ppos=ipad, h=h2, title=pad_title[ipad], 
                               xaxis_title='X [cm]', yaxis_title='Y [cm]',
                               draw_options='0 colz')

            #title and save
            ctitle = 'Layer #'+str(layers[isd])+'       Wafer: U='+uvN[0]+', V='+uvN[1] 
            t = TPaveLabel(0.02, 0.94, 0.98, 0.99, ctitle)
            plot.setCanvasTitle(cpos=0, t=t)
            picname = directory +'/hex_'+uvN[0]+'_'+uvN[1]
            plot.save(cpos=0, name=picname)
