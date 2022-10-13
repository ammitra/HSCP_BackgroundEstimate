import ROOT
import numpy as np

def makeRegion(name, iHist):
    '''
    takes in a 20x20 HSPC ProbQvsIas histogram and returns the 20x1 Pass and Fail histos
    '''
    # first do Pass
    hPassTemp = ROOT.TH2F('hpass','ProbQvsIas;ProbQ;Ias',1,0.0,1.0,20,0.0,1.0) # need to change x [0,0.1]?
    # loop over 
    for yBin in range(1,21):
        ySum = 0.0
        for xBin in range(1,3):
            print('Getting bin content: ({},{})'.format(xBin,yBin))
            ySum += iHist.GetBinContent(xBin,yBin)
        hPassTemp.SetBinContent(1,yBin,ySum)

    # now do Fail
    hFailTemp = ROOT.TH2F('hfail','ProbQvsIas;ProbQ;Ias',1,0.0,1.0,20,0.0,1.0) # need to change x [0.1,0.7]?
    for yBin in range(1,21):
        ySum = 0.0
        for xBin in range(3,15):
            print('Getting bin content: ({},{})'.format(xBin,yBin))
            ySum += iHist.GetBinContent(xBin,yBin)
        hFailTemp.SetBinContent(1,yBin,ySum)


    f = ROOT.TFile.Open('HSPC_{}.root'.format(name),'RECREATE')
    # 2DAlphabet expects signal/blinding to be on X axis only.. So need to invert
    hPass = ROOT.TH2F('hpass','IasVsProbQ;Ias;ProbQ',20,0.0,1.0,1,0.0,1.0) # need to change y [0.0, 0.1]?
    hFail = ROOT.TH2F('hfail','IasVsProbQ;Ias;ProbQ',20,0.0,1.0,1,0.0,1.0) # need to change y [0.1, 0.7]?
    for yBin in range(1,21):
	hPNew = hPassTemp.GetBinContent(1,yBin)
	hPass.SetBinContent(yBin,1,hPNew)
	hFNew = hFailTemp.GetBinContent(1,yBin)
	hFail.SetBinContent(yBin,1,hFNew)

    hFail.Write()
    hPass.Write()
    f.Close()

#f = ROOT.TFile.Open('crab_Analysis_2018_HSCPgluino_M-1800_CodeV40p4_v1.root')
#h = f.Get('HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas')

d = {'Signal': 'crab_Analysis_2018_HSCPgluino_M-1800_CodeV40p4_v1.root',
 'Data': 'crab_Analysis_2018_AllBackground_CodeV40p4_v1.root'}

hName = 'HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas'

for name, fName in d.items():
    fTemp = ROOT.TFile.Open(fName)
    h = fTemp.Get(hName)
    makeRegion(name, h)
    fTemp.Close()


