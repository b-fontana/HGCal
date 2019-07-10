import ROOT
from array import array

def buildMedianProfile(h):
    """build the median profile"""
    medianGr = ROOT.TGraphAsymmErrors()
    medianGr.SetName('{}_medianprof'.format(h.GetName()))
    medianGr.SetLineWidth(2)
    medianGr.SetMarkerStyle(20)

    #median and 1 sigma quantiles
    xq = array('d', [0.16,0.5,0.84]) 
    #storage for the x values that correspond to the above quantiles
    yq = array('d', [0.0 ,0.0,0.0 ])

    for xbin in xrange(1,h.GetNbinsX()+1):
        tmp=h.ProjectionY('tmp',xbin,xbin)
        tmp.GetQuantiles(3,yq,xq)
        xcen=h.GetXaxis().GetBinCenter(xbin)

        npts=medianGr.GetN()
        medianGr.SetPoint(npts,xcen,yq[1])
        medianGr.SetPointError(npts,0,0,yq[1]-yq[0],yq[2]-yq[1])
        tmp.Delete()
    return medianGr
