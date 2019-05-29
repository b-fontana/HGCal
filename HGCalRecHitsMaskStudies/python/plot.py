import sys
import os

import ROOT
import numpy as np
from matplotlib import pyplot as plt
import argparse

import RootUtils
import SystemUtils
import Argparser

def getUVFromTitle(title):
    u, v, rest = title.split(',')
    return u, v, rest

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
if len(histo_names_rechits) != len(histo_names_geom):
    raise ValueError('The files with the wafer codes need to have the same size!')
for _l in range(len(layers)):
    if len(histo_names_rechits[_l]) != len(histo_names_geom[_l]):
        raise ValueError('The files with the wafer codes need to have the same size!')

myfile = ROOT.TFile.Open(FLAGS.root_file+".root", "READ")
online_dir = '/eos/user/b/bfontana/www/'
pcoords = [[[0.02,0.49,0.98,0.93],    #pad1
           [0.02,0.02,0.98,0.48]]]   #pad2
for isd,sd in enumerate(subd_names):
    directory = os.path.join(online_dir,FLAGS.directory,sd)
    SystemUtils.createDir(directory)
    for nr, ng in zip(histo_names_rechits[isd],histo_names_geom[isd]):
        with RootUtils.RootPlotting(ncanvas=1, npads=2, pcoords=pcoords) as plot:
            uv1 = getUVFromTitle(nr)
            uv2 = getUVFromTitle(ng)
            if uv1[0] != uv2[0] or uv1[1] != uv2[1]:
                raise ValueError('The stored names have to match!')
            title = 'Layer #'+str(layers[isd])+'       Wafer: U='+uv1[0]+', V='+uv1[1] 
            print(os.path.join(maind,sd,nr))

            h1 = myfile.Get(os.path.join(maind,sd,nr))
            plot.plotHistogram(cpos=0, ppos=0, h=h1, 
                               title='RecHits', xaxis_title='u', yaxis_title='v')
            h2 = myfile.Get(os.path.join(maind,sd,ng))
            plot.plotHistogram(cpos=0, ppos=1, h=h2, 
                               title='Geometry', xaxis_title='u', yaxis_title='v')
            ttt = ROOT.TPaveLabel(0.02, 0.94, 0.98, 0.99, title)
            plot.setCanvasTitle(cpos=0, t=ttt)
            name = directory +'/pic_'+uv1[0]+'_'+uv1[1]
            plot.save(cpos=0, name=name)
