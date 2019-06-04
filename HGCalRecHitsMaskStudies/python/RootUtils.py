from ROOT import TCanvas, TPad, TPaveLabel, TPaveText, gSystem, gDirectory, kBlue
import numpy as np

class RootPlotting:
    def __init__(self, ncanvas, npads=[1], pcoords=None):
        """
        Arguments: 
        -> ncanvas: number of canvases that the instance will have
        -> npads: list number of pads per canvas, where each element in the list refers to a 
           different canvas
        -> pcoords: list of list of the coordinates of the pads, where each sublist refers to a different canvas
        """
        if ncanvas<=0: 
            raise ValueError('The object needs to possess at least one canvas')
        if type(npads) is int:
            npads = [npads]
        if any(_i<=0 for _i in npads): 
            raise ValueError('The object needs to possess at least one pad')        
        if any(_i==1 for _i in npads) and pcoords is None:
            raise ValueError('Please specify the coordinates associated to each pad.')
        if pcoords is not None and len(pcoords) != ncanvas:
            raise ValueError('The dimensions of your coordinate vector do not match the'
                             'number of canvases that was specified.')
        for _ic in range(ncanvas):
            if len(pcoords[_ic]) != npads[_ic]:
                raise ValueError('')

        self._c = [TCanvas('c','c',200,10,700,900) for _ in range(ncanvas)]
        self._p = [[] for _ in range(ncanvas)]
        if pcoords == None:
            for _ic in range(ncanvas):
                for _ip in range(npads[_ic]):
                    self._c[_ic].cd()
                    self._p[_ic].append(TPad('pad'+str(_ic)+'_'+str(_ip),'pad'+str(_ic)+'_'+str(_ip),
                                             0.03,0.03,0.97,0.8,21))
        else:
            for _ic in range(ncanvas):
                for _ip in range(npads[_ic]):
                    self._c[_ic].cd()
                    self._p[_ic].append(TPad('pad'+str(_ic)+'_'+str(_ip),'pad'+str(_ic)+'_'+str(_ip),
                                             pcoords[_ic][_ip][0],pcoords[_ic][_ip][1],
                                             pcoords[_ic][_ip][2],pcoords[_ic][_ip][3],21))
        for _ic in range(ncanvas):
            for _ip in range(npads[_ic]):
                self._c[_ic].cd()
                self._p[_ic][_ip].Draw()
                self._p[_ic][_ip].cd()
                self._p[_ic][_ip].GetFrame().SetFillColor( 18 )

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
                                   coords[0],coords[1],coords[2],coords[3],21)
            pad[cpos][ppos].GetFrame().SetFillColor( 18 )
            pad[cpos][ppos].Draw()
        except IndexError():
            print('The specified pad index does not correspond to an existing pad in canvas'
                  'number {}'.format(cpos))
            raise
    
    def setCanvasTitle(self, cpos, t):
        self._c[cpos].cd()
        t.SetFillColor(16)
        t.SetTextFont(52)
        t.Draw()

    def plotHistogram(self, cpos, ppos, h, 
                      title='', xaxis_title='', yaxis_title=''):
        """
        Arguments:
        -> cpos: canvas index according to constructor or later addition
        -> ppos: pad index inside a canvas according to constructor or later addition
        -> h: Root histogram to plot
        """    
        self._c[cpos].cd()
        self._p[cpos][ppos].cd()        
        h.SetTitle(title)
        h.GetXaxis().SetTitle(xaxis_title)
        h.GetYaxis().SetTitle(yaxis_title)
        h.Draw('colz')


    def save(self, cpos, name):
        self._c[cpos].SaveAs(name+'.png')

    def saveAll(self, name):
        for _i in range(len(self._c)):
            self.save(_i, name+"_v"+str(_i))
