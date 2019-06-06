import sys
import os

import ROOT
import numpy as np
from matplotlib import pyplot as plt
import argparse
import math, copy

import RootUtils
import SystemUtils
import Argparser

def getUVFromTitle(title):
    u, v, N, rest = title.split(',')
    return u, v, N, rest

def convertToHex(h, N):
    """
    Performs the conversion between standard (U,V) TH2D to
    (X,Y) hexagonal histogram within an hexagonal grid.
    """
    def addBins(h, xstart, ystart, a, k, s):
        from array import array
        x = array('d', [0., 0., 0., 0., 0., 0.])
        y = array('d', [0., 0., 0., 0., 0., 0.])
        xloop = xstart
        yloop = ystart + a*math.sqrt(3)/2.0;

        for sCounter in range(s):
            print("S beginning: ", sCounter)
            xtemp = xloop # Resets the temp variable
            
            #Determine the number of hexagons in that row
            if sCounter%2 == 0:
                numberOfHexagonsInTheRow = k
            else:
                numberOfHexagonsInTheRow = k - 1
 
            for kCounter in range(numberOfHexagonsInTheRow):
                print(sCounter, kCounter)
                #Go around the hexagon
                x[0] = xtemp
                y[0] = yloop
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

                h.AddBin(6, x, y)
 
                # Go right
                xtemp += 3*a
 
            # Increment the starting position
            if sCounter%2 == 0:
                xloop += 1.5*a
            else:
                xloop -= 1.5*a
            yloop += a*math.sqrt(3)/2
            print("S end: ", sCounter)            

    def convertCoords(u, v):
        return 1.5*(v-N)+1, u-0.5*(N+v)
    
    N2 = 2*N+1

    x00, y00 = convertCoords(0, 0)
    x01, y01 = convertCoords(0, 1)
    dist_side = N * math.sqrt( (x01-x00)**2 + (y01-y00)**2 )
    #bottom left graph corner
    startx = - dist_side*math.sin(math.pi/3)
    starty = y00 - dist_side*math.cos(math.pi/3)

    _, miny = convertCoords(0, N-1)
    _, maxy = convertCoords(N2, N-1)
    maxx, _ = convertCoords(N2, N2)

    hhex = ROOT.TH2Poly('hhex', 'hhex', startx-5, maxx, miny-5, maxy+15)
    addBins(hhex,startx,starty,(x01-x00)/math.sqrt(3)/100,N*100,2*N2*100)
    print(N2)
    """
    hhex = ROOT.TH2Poly('hhex', 'hhex', startx-5, startx+5, miny-5, miny+5)
    hhex.SetTitle("Honeycomb example")
    hhex.Honeycomb(startx,starty,(x01-x00)/2.,N2,N2)
    """

    for u in range(h.GetNbinsX()):
        for v in range(h.GetNbinsY()):
            x, y = convertCoords(u, v)   
            global_bin = hhex.FindBin(x, y)
            hhex.SetBinContent(global_bin, h.GetBinContent(u, v))

    return hhex

parser = argparse.ArgumentParser()
FLAGS, _ = Argparser.add_args(parser)
Argparser.print_args(FLAGS)

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(ROOT.kFALSE)
maind = 'outpath'
layers = [x for x in range(1,FLAGS.nlayers+1)]
subd_names = ['layer'+str(layers[x]) for x in range(len(layers))]
histo_names_rechits = [[] for _ in range(len(layers))]
histo_names_geom = [[] for _ in range(len(layers))]
for il,l in enumerate(layers):
    with open('data/HistoNamesRecHitsLayer'+str(l)+'_mask'+str(FLAGS.mask)+'.txt') as f:
        for line in f:
            histo_names_rechits[il].append(line[:-1])
    with open('data/HistoNamesGeomLayer'+str(l)+'_mask'+str(FLAGS.mask)+'.txt') as f:
        for line in f:
            histo_names_geom[il].append(line[:-1])

myfile = ROOT.TFile.Open(FLAGS.root_file, "READ")
online_dir = '/eos/user/b/bfontana/www/'
pcoords = [[[0.02,0.49,0.98,0.93],    #pad1
           [0.02,0.02,0.98,0.48]]]   #pad2
for isd,sd in enumerate(subd_names):
    directory = os.path.join(online_dir,FLAGS.directory,sd)
    SystemUtils.createDir(os.path.join(online_dir,FLAGS.directory))
    SystemUtils.createDir(directory)
    for nr,ng in zip(histo_names_rechits[isd],histo_names_geom[isd]):
        with RootUtils.RootPlotting(ncanvas=1, npads=2, pcoords=pcoords) as plot:
            #first pad
            uvN1 = getUVFromTitle(nr)
            h1 = myfile.Get(os.path.join(maind,sd,nr))
            h1 = convertToHex(h1, int(uvN1[2]))
            plot.plotHistogram(cpos=0, ppos=0, h=h1, 
                               title='RecHits', xaxis_title='X', yaxis_title='Y')
           
            #second pad
            uvN2 = getUVFromTitle(ng)
            h2 = myfile.Get(os.path.join(maind,sd,ng))
            h2 = convertToHex(h2, int(uvN2[2]))
            plot.plotHistogram(cpos=0, ppos=1, h=h2, 
                               title='Geometry', xaxis_title='X', yaxis_title='Y')

            #title and save
            title = 'Layer #'+str(layers[isd])+'       Wafer: U='+uvN1[0]+', V='+uvN1[1] 
            t = ROOT.TPaveLabel(0.02, 0.94, 0.98, 0.99, title)
            plot.setCanvasTitle(cpos=0, t=t)
            name = directory +'/hex_'+uvN1[0]+'_'+uvN1[1]
            plot.save(cpos=0, name=name)
            quit()
