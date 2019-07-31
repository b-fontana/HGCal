import abc
from array import array as Carray
from ROOT import TFile, TCanvas, TPad, TMath, TStyle, TExec, TLatex, TGraph, TF1
from ROOT import gSystem, gDirectory, gStyle, kBlue

def standardPadCoords(rows):
    if rows==2:
        return [[[0.01,0.51,0.33,0.99],
                 [0.34,0.51,0.66,0.99],  
                 [0.67,0.51,0.99,0.99],
                 [0.01,0.01,0.33,0.49],   
                 [0.34,0.01,0.66,0.49],   
                 [0.67,0.01,0.99,0.49]]]
    elif rows==4:
        return [[[0.01,0.755,0.33,0.995],   
                 [0.34,0.755,0.66,0.995],   
                 [0.67,0.755,0.99,0.995],
                 [0.01,0.505,0.33,0.745],   
                 [0.34,0.505,0.66,0.745],
                 [0.67,0.505,0.99,0.745],
                 [0.01,0.255,0.33,0.495],
                 [0.34,0.255,0.66,0.495],
                 [0.67,0.255,0.99,0.495],
                 [0.01,0.005,0.33,0.245],
                 [0.34,0.005,0.66,0.245],
                 [0.67,0.005,0.99,0.245]]]
    elif rows==5:
        return [[[0.01,0.805,0.33,0.995],   
                 [0.34,0.805,0.66,0.995],   
                 [0.67,0.805,0.99,0.995],
                 [0.01,0.605,0.33,0.795],   
                 [0.34,0.605,0.66,0.795],
                 [0.67,0.605,0.99,0.795],
                 [0.01,0.405,0.33,0.595],
                 [0.34,0.405,0.66,0.595],
                 [0.67,0.406,0.99,0.595],
                 [0.01,0.205,0.33,0.395],
                 [0.34,0.205,0.66,0.395],
                 [0.67,0.205,0.99,0.395],
                 [0.01,0.005,0.33,0.195],
                 [0.34,0.005,0.66,0.195],
                 [0.67,0.005,0.99,0.195]]]
    elif rows==6:
        return [[[0.01,0.835,0.33,0.995],   
                 [0.34,0.835,0.66,0.995],   
                 [0.67,0.835,0.99,0.995],
                 [0.01,0.670,0.33,0.830],   
                 [0.34,0.670,0.66,0.830],
                 [0.67,0.670,0.99,0.830],
                 [0.01,0.505,0.33,0.665],
                 [0.34,0.505,0.66,0.665],
                 [0.67,0.505,0.99,0.665],
                 [0.01,0.340,0.33,0.500],
                 [0.34,0.340,0.66,0.500],
                 [0.67,0.340,0.99,0.500],
                 [0.01,0.175,0.33,0.335],
                 [0.34,0.175,0.66,0.335],
                 [0.67,0.175,0.99,0.335],
                 [0.01,0.010,0.33,0.170],
                 [0.34,0.010,0.66,0.170],
                 [0.67,0.010,0.99,0.170]]]
    else:
        raise ValueError('standardPadCoords was not implemented for {} rows.'.format(rows))

class RootPlotting:
    def __init__(self, ncanvas, npads=[1], cdims=None, pcoords=None):
        """
        Manages Root plots.
        One canvas corresponds to one picture.
        The numbering of canvases and the number of pads per canvas starts at zero.
        Works for Python2.
        
        Arguments: 
        -> ncanvas: number of canvases that the instance will have
        -> npads: list number of pads per canvas, where each element in the list 
        refers to a different canvas
        -> cdims: list of lists of horizontal and vertical dimensions, where each
        sublist refers to a different canvas
        -> pcoords: list of list of the coordinates of the pads, where each sublist 
        refers to a different canvas
        """
        if ncanvas<=0: 
            raise ValueError('The object needs to possess at least one canvas')
        if type(npads) is int:
            npads = [npads]
        if any(_i<=0 for _i in npads): 
            raise ValueError('The object needs to possess at least one pad')        
        if any(_i>=1 for _i in npads) and pcoords is None:
            raise ValueError('Please specify the coordinates associated to each pad.')
        if pcoords is not None and len(pcoords) != ncanvas:
            raise ValueError('The dimensions of your coordinate vector do not match the '
                             'number of canvases that was specified.')
        if cdims is not None and len(cdims) != ncanvas:
            raise ValueError('The dimensions of your cdims vector do not match the '
                             'number of canvases that was specified.')
        
        self._c = []
        if cdims == None:
            for _ic in range(ncanvas):
                self._c.append(TCanvas('c'+str(_ic), 'c'+str(_ic), 200, 10, 
                                       500, 1000))
        else:
            for _ic in range(ncanvas):
                self._c.append(TCanvas('c'+str(_ic), 'c'+str(_ic), 200, 10, 
                                       cdims[_ic][0], cdims[_ic][1]))

        self._p = [[] for _ in range(ncanvas)]
        if pcoords == None:
            for _ic in range(ncanvas):
                for _ip in range(npads[_ic]):
                    self._c[_ic].cd()
                    self._p[_ic].append(TPad('pad'+str(_ic)+'_'+str(_ip),
                                             'pad'+str(_ic)+'_'+str(_ip),
                                             0.03,0.03,0.97,0.8,22))
        else:
            for _ic in range(ncanvas):
                for _ip in range(npads[_ic]):
                    self._c[_ic].cd()
                    self._p[_ic].append(TPad('pad'+str(_ic)+'_'+str(_ip),
                                             'pad'+str(_ic)+'_'+str(_ip),
                                             pcoords[_ic][_ip][0],pcoords[_ic][_ip][1],
                                             pcoords[_ic][_ip][2],pcoords[_ic][_ip][3],22))

        if len(cdims) != ncanvas:
            raise ValueError("The number of canvases is incompatible with the "
                             "the size of the 'cdims' vector")
        for _ic in range(ncanvas):
            if len(cdims[_ic]) != 2:
                raise ValueError('Two dimensions for each canvas have to be specified.')
            if len(pcoords[_ic]) != npads[_ic]:
                raise ValueError("The number of pads [{}] in canvas #{} is incompatible "
                                 "with the the size of the 'pcoords[{}]' [{}] entry."
                                 .format(len(pcoords[_ic]), _ic, _ic, npads[_ic]))

        for _ic in range(ncanvas):
            for _ip in range(npads[_ic]):
                self._c[_ic].cd()
                self._p[_ic][_ip].Draw()
                self._p[_ic][_ip].cd()
                self._p[_ic][_ip].GetFrame().SetFillColor( 16 )
                self._p[_ic][_ip].SetGrid()
                self._p[_ic][_ip].SetRightMargin(0.13);
                self._p[_ic][_ip].Draw()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        for _ic in range(len(self._c)):
            self._c[_ic].Close()
            gSystem.ProcessEvents()

    def setCanvas(self, cpos):
        try:
            self._c[cpos] = TCanvas('c','c',200,10,700,900)
        except IndexError:
            if cpos != len(self._c): raise
            else: self._c.append(TCanvas('c','c',200,10,700,900))
            
    def setPad(self, cpos, ppos, coords):
        try:
            self._c[cpos].cd()
        except IndexError():
            print('The specified canvas index does not correspond to an existing canvas')
            raise
        try:
            pad[cpos][ppos] = TPad('pad'+str(cpos)+'_'+str(ppos),
                                   'pad'+str(cpos)+'_'+str(ppos),
                                   coords[0],coords[1],coords[2],coords[3],22)
            pad[cpos][ppos].GetFrame().SetFillColor( 16 )
            pad[cpos][ppos].Draw()
        except IndexError():
            print('The specified pad index does not correspond to an existing pad in canvas'
                  'number {}'.format(cpos))
            raise
    
    def setCanvasTitle(self, cpos, t):
        self._c[cpos].cd()
        t.SetFillColor(16)
        t.SetTextFont(62)
        t.SetTextSize(.6)
        t.Draw()

    def setLatex(self, tf=42, ts=0.04):
        tex = TLatex()
        tex.SetTextFont(tf)
        tex.SetTextSize(ts)
        tex.SetNDC()
        return tex

    def setHistogram(self, cpos, ppos, h, 
                     title='', xaxis_title='', yaxis_title='', 
                     xranges=(None,None), yranges=(None,None),
                     lc=1, lw=3, mc=1, msize=1, mstyle=8):
        self._c[cpos].cd()
        self._p[cpos][ppos].cd()       
        h.SetTitle(title)
        h.SetTitleOffset(.5)
        if xaxis_title != '':
            h.GetXaxis().SetTitle(xaxis_title)
        h.GetXaxis().SetTitleOffset(1.2)
        if yaxis_title != '':
            h.GetYaxis().SetTitle(yaxis_title)
        h.GetYaxis().SetTitleOffset(1.2)
        if xranges != (None,None):
            h.GetXaxis().SetRangeUser(xranges[0],xranges[1]);
        if yranges != (None,None):
            h.GetYaxis().SetRangeUser(yranges[0],yranges[1]);
        h.SetLineColor(lc)
        h.SetMarkerColor(mc)
        h.SetMarkerSize(msize)
        h.SetMarkerStyle(mstyle)
        h.SetLineWidth(lw)


    def plotHistogram(self, cpos, ppos, h, 
                      draw_options='colz', copy=False, hset=True, **kwargs):
        """
        Arguments:
        -> cpos: canvas index according to constructor or later addition
        -> ppos: pad index inside a canvas according to constructor or later addition
        -> h: Root histogram to plot
        """
        if hset:
            self.setHistogram(cpos, ppos, h, **kwargs)
        if copy:
            h.DrawCopy(draw_options)
        else:
            h.Draw(draw_options)

    def fitHistogram(self, h, f=None, fname='crystalball', frange=(-5.,5.), tex=None):
        """
        Fits TF1 'f' to histogram 'h'. 
        If f is 'None', it uses one of the provided functions:
        -> crystalball
        -> gaussian
        """
        if f == None:
            f = TF1("f", fname, frange[0], frange[1])
            if fname == 'crystalball':
                f.SetParameters(1000., 0., .2, 2., 1.)
                f.SetParLimits(0, 1., 10000)
                f.SetParLimits(1, -.3, .3)
                f.SetParLimits(2, 0.00001, 1.)
                f.SetParLimits(3, 0.00001, 5.)
                f.SetParLimits(4, 0.00001, 5.)
            elif fname == 'gaus':
                f.SetParameters(1., 0., 1.)
            else:
                raise ValueError('The specified function is not supported.')
        h.Fit(f.GetName(), 'M+')
        if tex != 'None':
            tex.DrawLatex(0.15,0.84, '#mu={}#pm{}'
                          .format(round(f.GetParameter(1),4), round(f.GetParError(1),4)))
            tex.DrawLatex(0.15,0.8, '#sigma={}#pm{}'
                          .format(round(f.GetParameter(2),4), round(f.GetParError(2),4)))

    def setGraph(self, cpos, ppos, g,
                 title='', name='', xaxis_title='', yaxis_title='',
                 xranges=(None,None), yranges=(None,None),
                 lc=1, lw=3, mc=1, msize=1, mstyle=8,
                 draw_options='', copy=False):
        self._c[cpos].cd()
        self._p[cpos][ppos].cd()       
        g.SetTitle(title)
        g.SetName(name)
        if xaxis_title != '':
            g.GetXaxis().SetTitle(xaxis_title)
        if yaxis_title != '':
            g.GetYaxis().SetTitle(yaxis_title)
        if xranges != (None,None):
            g.GetXaxis().SetRangeUser(xranges[0],xranges[1]);
        if yranges != (None,None):
            g.GetYaxis().SetRangeUser(yranges[0],yranges[1]);
        g.GetXaxis().SetTitleOffset(1.2)
        g.GetYaxis().SetTitleOffset(1.2)
        g.SetLineColor(lc)
        g.SetMarkerColor(mc)
        g.SetMarkerSize(msize)
        g.SetMarkerStyle(mstyle)
        g.SetLineWidth(lw)

    def plotGraph(self, cpos, ppos, g,
                  draw_options='AP', copy=False, gset=True, **kwargs): 
        if gset:
            self.setGraph(cpos, ppos, g, **kwargs)
        if copy:
            g.DrawCopy(draw_options)
        else:
            g.Draw(draw_options)

    def save(self, cpos, name):
        self._c[cpos].SaveAs(name)

    def saveAll(self, name):
        for _i in range(len(self._c)):
            self.save(_i, name+"_v"+str(_i))
