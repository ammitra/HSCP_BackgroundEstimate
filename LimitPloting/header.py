import ROOT
from ROOT import *
from contextlib import contextmanager
import os, pickle, subprocess, time,random
import math, collections
from math import sqrt
import array
import json
import CMS_lumi, tdrstyle
import pprint
pp = pprint.PrettyPrinter(indent = 2)

def setSnapshot(d=''):
    # header.executeCmd('combine -M MultiDimFit -d '+base_workspace+' --saveWorkspace --freezeParameters r --setParameters r=0,'+mask_string)
    # f = TFile.Open('higgsCombineTest.MultiDimFit.mH120.root')
    # w = f.Get('w')
    # w.loadSnapshot("MultiDimFit")
    # w.var("r").setConstant(0)
    # w.var("r").setVal(0)
    # w.var("r").setMin(-20)
    # w.var("r").setMax(20)
    # myargs = RooArgSet(w.allVars())
    # myargs.add(w.allCats())
    # w.saveSnapshot("initialFit",myargs)
    # fout = TFile('initialFitWorkspace.root',"recreate")
    # fout.WriteTObject(w,'w')
    # fout.Close()
    w_f = TFile.Open(d+'higgsCombineTest.FitDiagnostics.mH120.root')
    w = w_f.Get('w')
    fr_f = TFile.Open(d+'fitDiagnostics.root')
    fr = fr_f.Get('fit_b')
    myargs = RooArgSet(fr.floatParsFinal())
    w.saveSnapshot('initialFit',myargs,True)
    fout = TFile('initialFitWorkspace.root',"recreate")
    fout.WriteTObject(w,'w')
    fout.Close()

# Function stolen from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii
def openJSON(f,twoDconfig=True):
    with open(f) as fInput_config:
        input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

        if twoDconfig:
            for process in [proc for proc in input_config['PROCESS'].keys() if proc != 'HELP']:
                for index,item in enumerate(input_config['PROCESS'][process]['SYSTEMATICS']):           # There's one list that also
                    input_config['PROCESS'][process]['SYSTEMATICS'][index] = item.encode('ascii')  

    return input_config

def ascii_encode_dict(data):    
    ascii_encode = lambda x: x.encode('ascii') if isinstance(x, unicode) else x 
    return dict(map(ascii_encode, pair) for pair in data.items())

def copyHistWithNewXbins(thisHist,newXbins,copyName):
    # Make a copy with the same Y bins but new X bins
    ybins = []
    for iy in range(1,thisHist.GetNbinsY()+1):
        ybins.append(thisHist.GetYaxis().GetBinLowEdge(iy))
    ybins.append(thisHist.GetYaxis().GetXmax())

    ybins_array = array.array('f',ybins)
    ynbins = len(ybins_array)-1

    xbins_array = array.array('f',newXbins)
    xnbins = len(xbins_array)-1 

    # Use copyName with _temp to avoid overwriting if thisHist has the same name
    # We can do this at the end but not before we're finished with thisHist
    hist_copy = TH2F(copyName+'_temp',copyName+'_temp',xnbins,xbins_array,ynbins,ybins_array)
    hist_copy.Sumw2()
    
    hist_copy.GetXaxis().SetName(thisHist.GetXaxis().GetName())
    hist_copy.GetYaxis().SetName(thisHist.GetYaxis().GetName())

    # Loop through the old bins
    for ybin in range(1,ynbins+1):
        # print 'Bin y: ' + str(binY)
        for xbin in range(1,xnbins+1):
            new_bin_content = 0
            new_bin_errorsq = 0
            new_bin_min = hist_copy.GetXaxis().GetBinLowEdge(xbin)
            new_bin_max = hist_copy.GetXaxis().GetBinUpEdge(xbin)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for old_xbin in range(1,thisHist.GetNbinsX()+1):
                old_bin_min = thisHist.GetXaxis().GetBinLowEdge(old_xbin)
                old_bin_max = thisHist.GetXaxis().GetBinUpEdge(old_xbin)
                if old_bin_min >= new_bin_max:
                    break
                elif old_bin_min >= new_bin_min and old_bin_min < new_bin_max:
                    if old_bin_max <= new_bin_max:
                        new_bin_content += thisHist.GetBinContent(old_xbin,ybin)
                        new_bin_errorsq += thisHist.GetBinError(old_xbin,ybin)**2
                    elif old_bin_max > new_bin_max:
                        raise ValueError(
                            '''The requested X rebinning does not align bin edges with the input bin edge.
                            Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(old_bin_min,old_bin_max,new_bin_min,new_bin_max))
                elif old_bin_min <= new_bin_min and old_bin_max > new_bin_min:
                    raise ValueError(
                        '''The requested Y rebinning does not align bin edges with the input bin edge.
                        Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(old_bin_min,old_bin_max,new_bin_min,new_bin_max))

            # print '\t Setting content ' + str(newBinContent) + '+/-' + str(sqrt(newBinErrorSq))
            if new_bin_content > 0:
                hist_copy.SetBinContent(xbin,ybin,new_bin_content)
                hist_copy.SetBinError(xbin,ybin,sqrt(new_bin_errorsq))

    # Will now set the copyName which will overwrite thisHist if it has the same name
    hist_copy.SetName(copyName)
    hist_copy.SetTitle(copyName)

    return hist_copy

def copyHistWithNewYbins(thisHist,newYbins,copyName):
    # Make a copy with the same X bins but new Y bins
    xbins = []
    for ix in range(1,thisHist.GetNbinsX()+1):
        xbins.append(thisHist.GetXaxis().GetBinLowEdge(ix))
    xbins.append(thisHist.GetXaxis().GetXmax())

    xbins_array = array.array('f',xbins)
    xnbins = len(xbins_array)-1

    ybins_array = array.array('f',newYbins)
    ynbins = len(ybins_array)-1

    # Use copyName with _temp to avoid overwriting if thisHist has the same name
    # We can do this at the end but not before we're finished with thisHist
    hist_copy = TH2F(copyName+'_temp',copyName+'_temp',xnbins,xbins_array,ynbins,ybins_array)
    hist_copy.Sumw2()
    
    hist_copy.GetXaxis().SetName(thisHist.GetXaxis().GetName())
    hist_copy.GetYaxis().SetName(thisHist.GetYaxis().GetName())

    # Loop through the old bins
    for xbin in range(1,xnbins+1):
        # print 'Bin y: ' + str(binY)
        for ybin in range(1,ynbins+1):
            new_bin_content = 0
            new_bin_errorsq = 0
            new_bin_min = hist_copy.GetYaxis().GetBinLowEdge(ybin)
            new_bin_max = hist_copy.GetYaxis().GetBinUpEdge(ybin)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for old_ybin in range(1,thisHist.GetNbinsY()+1):
                old_bin_min = thisHist.GetYaxis().GetBinLowEdge(old_ybin)
                old_bin_max = thisHist.GetYaxis().GetBinUpEdge(old_ybin)
                if old_bin_min >= new_bin_max:
                    break
                elif old_bin_min >= new_bin_min and old_bin_min < new_bin_max:
                    if old_bin_max <= new_bin_max:
                        new_bin_content += thisHist.GetBinContent(xbin,old_ybin)
                        new_bin_errorsq += thisHist.GetBinError(xbin,old_ybin)**2
                    elif old_bin_max > new_bin_max:
                        raise ValueError(
                            '''The requested Y rebinning does not align bin edges with the input bin edge.
                            Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(old_bin_min,old_bin_max,new_bin_min,new_bin_max))
                elif old_bin_min <= new_bin_min and old_bin_max > new_bin_min:
                    raise ValueError(
                        '''The requested Y rebinning does not align bin edges with the input bin edge.
                        Cannot split input bin [%s,%s] with output bin [%s,%s]'''%(old_bin_min,old_bin_max,new_bin_min,new_bin_max))

            # print '\t Setting content ' + str(newBinContent) + '+/-' + str(sqrt(newBinErrorSq))
            if new_bin_content > 0:
                hist_copy.SetBinContent(xbin,ybin,new_bin_content)
                hist_copy.SetBinError(xbin,ybin,sqrt(new_bin_errorsq))

    # Will now set the copyName which will overwrite thisHist if it has the same name
    hist_copy.SetName(copyName)
    hist_copy.SetTitle(copyName)

    return hist_copy

def ConvertToEvtsPerUnit(hist,width=None):
    if width == None:
        use_width = GetMinWidth(hist)
    else:
        use_width = width

    converted = hist.Clone()
    for ibin in range(1,hist.GetNbinsX()+1):
        if hist.GetBinWidth(ibin) == use_width:
            continue
        else:
            factor = use_width/hist.GetBinWidth(ibin)
            new_content = factor * converted.GetBinContent(ibin)
            new_error = factor * converted.GetBinError(ibin)
            converted.SetBinContent(ibin,new_content)
            converted.SetBinError(ibin,new_error)
    
    return converted

def GetMinWidth(hist):
    use_width = 10**6
    for ibin in range(1,hist.GetNbinsX()+1):
        if hist.GetBinWidth(ibin) < use_width:
            use_width = hist.GetBinWidth(ibin)
    return int(use_width)

def smoothHist2D(name,histToSmooth,renormalize=False,iterate=1,skipEdges=False):
    print "Smoothing "+name
    if renormalize: norm = histToSmooth.Integral()
    smoothed_hist = histToSmooth.Clone(name)
    smoothed_hist.Reset()
    smoothed_hist.Sumw2()

    for ix in range(1,histToSmooth.GetNbinsX()+1):
        for iy in range(1,histToSmooth.GetNbinsY()+1):
            bin_contents = [histToSmooth.GetBinContent(ix,iy)]

            if ix == 1:
                if iy == 1: # lower left corner
                    if skipEdges: pass
                    else:
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                elif iy == histToSmooth.GetNbinsY(): # upper left corner
                    if skipEdges: pass
                    else:
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy-1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))
                else: # left wall
                    if not skipEdges:
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy-1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))
            
            elif ix == histToSmooth.GetNbinsX():
                if iy == 1: # lower right corner
                    if skipEdges: pass
                    else:
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                elif iy == histToSmooth.GetNbinsY(): # upper right corner
                    if skipEdges: pass
                    else:
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy-1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))
                else: # right wall
                    if not skipEdges: 
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy-1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))

            else:
                if iy == 1: # bottom wall
                    if not skipEdges:
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                    bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                elif iy == histToSmooth.GetNbinsY(): # top wall
                    if not skipEdges:
                        bin_contents.append(histToSmooth.GetBinContent(ix-1,iy-1))
                        bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))
                        bin_contents.append(histToSmooth.GetBinContent(ix+1,iy-1))
                    bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                    bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                else: # full square
                    bin_contents.append(histToSmooth.GetBinContent(ix+1,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix-1,iy+1))
                    bin_contents.append(histToSmooth.GetBinContent(ix+1,iy  ))
                    # bin_contents.append(histToSmooth.GetBinContent(ix  ,iy  ))
                    bin_contents.append(histToSmooth.GetBinContent(ix-1,iy  ))
                    bin_contents.append(histToSmooth.GetBinContent(ix+1,iy-1))
                    bin_contents.append(histToSmooth.GetBinContent(ix  ,iy-1))
                    bin_contents.append(histToSmooth.GetBinContent(ix-1,iy-1))

            avg = 0
            for b in bin_contents: avg += b
            avg = avg/len(bin_contents)

            smoothed_hist.SetBinContent(ix,iy,avg)

    if renormalize:
        smoothed_norm = smoothed_hist.Integral()
        smoothed_hist.Scale(norm/smoothed_norm)
    if iterate > 1: smoothed_hist = smoothHist2D(name+str(iterate),smoothed_hist,iterate=iterate-1)
    return smoothed_hist

def zeroNegativeBins(name,inhist):
    outhist = inhist.Clone(name)
    outhist.Reset()
    for ix in range(1,inhist.GetNbinsX()):
        for iy in range(1,inhist.GetNbinsY()):
            content = max(0,inhist.GetBinContent(ix,iy))
            outhist.SetBinContent(ix,iy,content)
            if content == 0: outhist.SetBinError(ix,iy,0)
            else: outhist.SetBinError(ix,iy,inhist.GetBinError(ix,iy))

    return outhist

def stitchHistsInX(name,xbins,ybins,thisHistList,blinded=[]):
    # Required that thisHistList be in order of desired stitching
    # `blinded` is a list of the index of regions you wish to skip/blind
    axbins = array.array('d',xbins)
    aybins = array.array('d',ybins)
    stitched_hist = TH2F(name,name,len(xbins)-1,axbins,len(ybins)-1,aybins)

    bin_jump = 0
    for i,h in enumerate(thisHistList):
        if i in blinded:
            bin_jump += thisHistList[i].GetNbinsX()
            continue
        
        for ybin in range(1,h.GetNbinsY()+1):
            for xbin in range(1,h.GetNbinsX()+1):
                stitched_xindex = xbin + bin_jump

                stitched_hist.SetBinContent(stitched_xindex,ybin,h.GetBinContent(xbin,ybin))
                stitched_hist.SetBinError(stitched_xindex,ybin,h.GetBinError(xbin,ybin))

        bin_jump += thisHistList[i].GetNbinsX()

    return stitched_hist

def rebinY(thisHist,name,tag,new_y_bins_array):
    xnbins = thisHist.GetNbinsX()
    xmin = thisHist.GetXaxis().GetXmin()
    xmax = thisHist.GetXaxis().GetXmax()

    rebinned = TH2F(name,name,xnbins,xmin,xmax,len(new_y_bins_array)-1,new_y_bins_array)
    # print new_y_bins_array
    # print rebinned.GetYaxis().GetBinUpEdge(len(new_y_bins_array)-1)
    rebinned.Sumw2()

    for xbin in range(1,xnbins+1):
        newBinContent = 0
        newBinErrorSq = 0
        rebinHistYBin = 1
        nybins = 0
        for ybin in range(1,thisHist.GetNbinsY()+1):
            # If upper edge of old Rpf ybin is < upper edge of rebinHistYBin then add the Rpf bin to the count
            if thisHist.GetYaxis().GetBinUpEdge(ybin) < rebinned.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += thisHist.GetBinContent(xbin,ybin)
                newBinErrorSq += thisHist.GetBinError(xbin,ybin)**2
                nybins+=1
            # If ==, add to newBinContent, assign newBinContent to current rebinHistYBin, move to the next rebinHistYBin, and restart newBinContent at 0
            elif thisHist.GetYaxis().GetBinUpEdge(ybin) == rebinned.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += thisHist.GetBinContent(xbin,ybin)
                newBinErrorSq += thisHist.GetBinError(xbin,ybin)**2
                nybins+=1
                rebinned.SetBinContent(xbin, rebinHistYBin, newBinContent/float(nybins))
                rebinned.SetBinError(xbin, rebinHistYBin, sqrt(newBinErrorSq)/float(nybins))# NEED TO SET BIN ERRORS
                rebinHistYBin += 1
                newBinContent = 0
                newBinErrorSq = 0
                nybins = 0
            else:
                print 'ERROR when doing psuedo-2D y rebin approximation. Slices do not line up on y bin edges'
                print 'Input bin upper edge = '+str(thisHist.GetYaxis().GetBinUpEdge(ybin))
                print 'Rebin upper edge = '+str(rebinned.GetYaxis().GetBinUpEdge(rebinHistYBin))
                quit()

    makeCan(name+'_rebin_compare',tag,[rebinned,thisHist])
    return rebinned

def splitBins(binList, sigLow, sigHigh):
    return_bins = {'LOW':[],'SIG':[],'HIGH':[]}
    for b in binList:
        if b <= sigLow:
            return_bins['LOW'].append(b)
        if b >= sigLow and b <= sigHigh:
            return_bins['SIG'].append(b)
        if b >= sigHigh:
            return_bins['HIGH'].append(b)

    return return_bins 

def remapToUnity(hist):
    # Map to [-1,1]
    ybins = array.array('d',[(hist.GetYaxis().GetBinLowEdge(b)-hist.GetYaxis().GetXmin())/(hist.GetYaxis().GetXmax()-hist.GetYaxis().GetXmin()) for b in range(1,hist.GetNbinsY()+1)]+[1])
    xbins = array.array('d',[(hist.GetXaxis().GetBinLowEdge(b)-hist.GetXaxis().GetXmin())/(hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin()) for b in range(1,hist.GetNbinsX()+1)]+[1])

    remap = TH2F(hist.GetName()+'_unit',hist.GetName()+'_unit',hist.GetNbinsX(),xbins,hist.GetNbinsY(),ybins)
    remap.Sumw2()

    for xbin in range(hist.GetNbinsX()+1):
        for ybin in range(hist.GetNbinsY()+1):
            remap.SetBinContent(xbin,ybin,hist.GetBinContent(xbin,ybin))
            remap.SetBinError(xbin,ybin,hist.GetBinError(xbin,ybin))

    return remap

def makeBlindedHist(nomHist,sigregion):
    # Grab stuff to make it easier to read
    xlow = nomHist.GetXaxis().GetXmin()
    xhigh = nomHist.GetXaxis().GetXmax()
    xnbins = nomHist.GetNbinsX()
    ylow = nomHist.GetYaxis().GetXmin()
    yhigh = nomHist.GetYaxis().GetXmax()
    ynbins = nomHist.GetNbinsY()
    blindName = nomHist.GetName()

    # Need to change nominal hist name or we'll get a memory leak
    nomHist.SetName(blindName+'_unblinded')

    blindedHist = TH2F(blindName,blindName,xnbins,xlow,xhigh,ynbins,ylow,yhigh)
    blindedHist.Sumw2()

    for binY in range(1,ynbins+1):
        # Fill only those bins outside the signal region
        for binX in range(1,xnbins+1):
            if nomHist.GetXaxis().GetBinUpEdge(binX) <= sigregion[0] or nomHist.GetXaxis().GetBinLowEdge(binX) >= sigregion[1]:
                if nomHist.GetBinContent(binX,binY) > 0:
                    blindedHist.SetBinContent(binX,binY,nomHist.GetBinContent(binX,binY))
                    blindedHist.SetBinError(binX,binY,nomHist.GetBinError(binX,binY))

    return blindedHist

def makeRDH(myTH2,RAL_vars,altname=''):
    name = myTH2.GetName()
    if altname != '':
        name = altname
    thisRDH = RooDataHist(name,name,RAL_vars,myTH2)
    return thisRDH

def makeRHP(myRDH,RAL_vars):
    name = myRDH.GetName()
    thisRAS = RooArgSet(RAL_vars)
    thisRHP = RooHistPdf(name,name,thisRAS,myRDH)
    return thisRHP

def colliMate(myString,width=18):
    sub_strings = myString.split(' ')
    new_string = ''
    for i,sub_string in enumerate(sub_strings):
        string_length = len(sub_string)
        n_spaces = width - string_length
        if i != len(sub_strings)-1:
            if n_spaces <= 0:
                n_spaces = 2
            new_string += sub_string + ' '*n_spaces
        else:
            new_string += sub_string
    return new_string

def dictStructureCopy(inDict):
    newDict = {}
    for k1,v1 in inDict.items():
        if type(v1) == dict:
            newDict[k1] = dictStructureCopy(v1)
        else:
            newDict[k1] = 0
    return newDict

def dictCopy(inDict):
    newDict = {}
    for k1,v1 in inDict.items():
        if type(v1) == dict:
            newDict[k1] = dictCopy(v1)
        else:
            newDict[k1] = v1
    return newDict

def printWorkspace(myfile,myworkspace):
    myf = TFile.Open(myfile)
    myw = myf.Get(myworkspace)
    myw.Print()

def ftestInfoLookup(projInfoDict):
    nrpfparams = 0
    nbins = 0
    for k in projInfoDict.keys():
        this_nrpfparams = len(projInfoDict[k]['rpfVarNames'])
        this_nbins = (len(projInfoDict[k]['full_x_bins'])-1) * (len(projInfoDict[k]['newYbins'])-1)
        if projInfoDict[k]['blindedFit'] == True:
            this_nbins = this_nbins - ( (len(projInfoDict[k]['newXbins']['SIG'])-1) * (len(projInfoDict[k]['newYbins'])-1) )

        nrpfparams+=this_nrpfparams
        nbins+=this_nbins
    
    return nrpfparams,nbins

def FStatCalc(filename1,filename2,p1,p2,n):
    print 'Calculating F statistic'
    # Flip flop to make sure p2 is always greater than p1 (more parameters should always fit better)
    if p1 > p2:
        p1, p2 = p2, p1
        filename1, filename2 = filename2, filename1
    print 'Files: ',filename1,filename2
    print 'Parameters: p1 %f, p2 %f, n %f'%(p1,p2,n)

    # Get limit trees from each file
    file1 = TFile.Open(filename1)
    tree1 = file1.Get("limit")
    file2 = TFile.Open(filename2)
    tree2 = file2.Get("limit")
    diffs=[]
    # print 'Entries ',tree1.GetEntries(),tree2.GetEntries()

    # Loop over entries and calculate the F statistics
    for i in range(0,tree1.GetEntries()):
        tree1.GetEntry(i)
        tree2.GetEntry(i)
        # print 'Limits ',tree1.limit,tree2.limit
        if tree1.limit-tree2.limit>0:
            F = (tree1.limit-tree2.limit)/(p2-p1)/(tree2.limit/(n-p2))
            # print 'Entry ',i, ":", tree2.limit, "-", tree1.limit, "=", tree2.limit-tree1.limit, "F =", F
            # if F < 50: 
            diffs.append(F)
        else:
            print 'WARNING in calculation of F statistic for entry %i. limit1-limit2 <=0 (%f - %f)' %(i,tree1.limit,tree2.limit)
            diffs.append(0)
    # print 'Diffs F stat: ',diffs
    return diffs

def makeToyCard(channelNames):
    # Open a new card
    toy_gen_card = open('card_toygen.txt','w')
    column_width = 20

    nchannels = len(channelNames) # imax
    nprocesses = 1 # jmax
    nsystematics = 0 # kmax

    toy_gen_card.write('imax '+str(nchannels)+'\n')      
    toy_gen_card.write('jmax '+str(nprocesses)+'\n')
    toy_gen_card.write('kmax '+str(nsystematics)+'\n')
    toy_gen_card.write('-'*120+'\n')

    processes = ['TotalSig','TotalBkg','data_obs']
    for proc in processes:
        if proc == 'TotalSig' or proc == 'data_obs':
            preORpost = '_prefit'
        else:
            preORpost = '_postfit'

        for chan in channelNames:
            toy_gen_card.write(colliMate('shapes  '+proc+' '+chan+' toy_gen_workspace.root w_toys:'+proc+'_'+chan+'\n'))

    toy_gen_card.write('-'*120+'\n')

    tempString = 'bin  '
    for chan in channelNames:
        tempString += (chan+' ')
    tempString += '\n'
    toy_gen_card.write(colliMate(tempString,column_width))

    tempString = 'observation  '
    for ichan in range(int(nchannels)):
        tempString += '-1 '
    tempString += '\n'
    toy_gen_card.write(colliMate(tempString,column_width))

    toy_gen_card.write('-'*120+'\n')

    bin_line = 'bin  '
    processName_line = 'process  '
    processCode_line = 'process  '
    rate_line = 'rate  '

    for chan in channelNames:
        for proc in ['TotalBkg','TotalSig']:
            # Start lines
            bin_line += (chan+' ')
            processName_line += (proc+' ')

            # If signal
            if proc == 'TotalSig':
                processCode_line += ('0 ')
                rate_line += ('-1 ')

            # If bkg
            elif proc == 'TotalBkg':
                processCode_line += ('1 ')
                rate_line += '-1 '                                            

    toy_gen_card.write(colliMate(bin_line+'\n',column_width))
    toy_gen_card.write(colliMate(processName_line+'\n',column_width))
    toy_gen_card.write(colliMate(processCode_line+'\n',column_width))
    toy_gen_card.write(colliMate(rate_line+'\n',column_width))
    toy_gen_card.write('-'*120+'\n')

    toy_gen_card.close()

def projInfoLookup(projDir,card_tag):
    # Check if there was more than one 2DAlphabet object run over
    more_than_one = False
    run_card = open(projDir+'/card_'+card_tag+'.txt','r')
    firstline = run_card.readline()
    if 'Combination of ' in firstline:
        more_than_one = True
    
    proj_info = {}
    if more_than_one:
        # Look up all categories run in the most recent fit
        twoD_names = getTwoDAlphaNames(firstline)
        for n in twoD_names:
            proj_info[n] = pickle.load(open(projDir+'/'+n+'/saveOut.p','r'))
            proj_info[n]['rpfVarNames'] = proj_info[n]['rpf'].getFuncVarNames()

    elif not more_than_one:
        proj_info[card_tag] = pickle.load(open(projDir+'/saveOut.p','r'))
        proj_info['rpfVarNames'] = proj_info['rpf'].getFuncVarNames()

    return proj_info

def getTwoDAlphaNames(line):
    card_locs = [loc for loc in line.split(' ') if (loc != ' ' and loc != '' and 'txt' in loc)]
    proj_names = [n.split('/')[0] for n in card_locs]
    return proj_names

def executeCmd(cmd,dryrun=False):
    print 'Executing: '+cmd
    if not dryrun:
        subprocess.call([cmd],shell=True)

def dictToLatexTable(dict2convert,outfilename,roworder=[],columnorder=[]):
    # First set of keys are row, second are column
    if len(roworder) == 0:
        rows = sorted(dict2convert.keys())
    else:
        rows = roworder
    if len(columnorder) == 0:
        columns = []
        for r in rows:
            thesecolumns = dict2convert[r].keys()
            for c in thesecolumns:
                if c not in columns:
                    columns.append(c)
        columns.sort()
    else:
        columns = columnorder

    latexout = open(outfilename,'w')
    latexout.write('\\begin{table}[] \n')
    latexout.write('\\begin{tabular}{|c|'+len(columns)*'c'+'|} \n')
    latexout.write('\\hline \n')

    column_string = ' &'
    for c in columns:
        column_string += str(c)+'\t& '
    column_string = column_string[:-2]+'\\\ \n'
    latexout.write(column_string)

    latexout.write('\\hline \n')
    for r in rows:
        row_string = '\t'+r+'\t& '
        for c in columns:
            if c in dict2convert[r].keys():
                row_string += str(dict2convert[r][c])+'\t& '
            else:
                row_string += '- \t& '
        row_string = row_string[:-2]+'\\\ \n'
        latexout.write(row_string)

    latexout.write('\\hline \n')
    latexout.write('\\end{tabular} \n')
    latexout.write('\\end{table}')
    latexout.close()

def reorderHists(histlist):
    if len(histlist) != 6:
        print Exception('reorderHists() only built to rearrange list of six hists from 2x3 to 3x2')
        return histlist

    outlist = []
    outlist.append(histlist[0])
    outlist.append(histlist[3])
    outlist.append(histlist[1])
    outlist.append(histlist[4])
    outlist.append(histlist[2])
    outlist.append(histlist[5])

    return outlist

def makeCan(name, tag, histlist, bkglist=[],totalBkg=None,signals=[],colors=[],
            titles=[],subtitles=[],sliceVar='X',dataName='Data',bkgNames=[],signalNames=[],logy=False,
            rootfile=False,xtitle='',ytitle='',ztitle='',dataOff=False,
            datastyle='pe',year=1, addSignals=True, extraText=''):  
    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]

    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1800
        height = 600
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        height = 1600
        width = 1200
        padx = 2
        pady = 3
        histlist = reorderHists(histlist)
        if bkglist != []: bkglist = reorderHists(bkglist)
        if signals != []: signals = reorderHists(signals)
        if totalBkg != None: totalBkg = reorderHists(totalBkg)
        if titles != []: titles = reorderHists(titles)
        if subtitles != []: subtitles = reorderHists(subtitles)
    else:
        print 'histlist of size ' + str(len(histlist)) + ' not currently supported'
        print histlist
        return 0

    tdrstyle.setTDRStyle()
    gStyle.SetLegendFont(42)
    gStyle.SetTitleBorderSize(0)
    gStyle.SetTitleAlign(33)
    gStyle.SetTitleX(.77)
        
    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [kRed,kMagenta,kGreen,kCyan,kBlue]
    if len(colors) == 0:   
        colors = default_colors
    color_idx_order = None
    stacks = []
    tot_hists_err = []
    tot_hists = []
    legends = []
    mains = []
    subs = []
    pulls = []
    logString = ''
    tot_sigs = []

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        # if len(histlist) > 1:
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd()        
        thisPad.SetRightMargin(0.0)
        thisPad.SetTopMargin(0.0)
        thisPad.SetBottomMargin(0.0)

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            gPad.SetLeftMargin(0.15)
            gPad.SetRightMargin(0.2)
            gPad.SetBottomMargin(0.12)
            gPad.SetTopMargin(0.1)
            if logy: gPad.SetLogz()
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            hist.GetZaxis().SetTitle(ztitle)
            hist.GetXaxis().SetTitleOffset(1.15)
            hist.GetYaxis().SetTitleOffset(1.5)
            hist.GetZaxis().SetTitleOffset(1.5)
            hist.GetYaxis().SetLabelSize(0.05)
            hist.GetYaxis().SetTitleSize(0.05)
            hist.GetXaxis().SetLabelSize(0.05)
            hist.GetXaxis().SetTitleSize(0.05)
            hist.GetZaxis().SetLabelSize(0.05)
            hist.GetZaxis().SetTitleSize(0.05)
            hist.GetXaxis().SetNdivisions(505)
            # hist.GetXaxis().SetLabelOffset(0.02)
            if 'lego' in datastyle.lower(): hist.GetZaxis().SetTitleOffset(1.4)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            if datastyle != 'pe': hist.Draw(datastyle)
            else: hist.Draw('colz')
            if len(bkglist) > 0:
                raise TypeError('ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.')

            CMS_lumi.extraText = extraText
            CMS_lumi.CMS_lumi(thisPad, year, 11, sim=False if 'data' in name.lower() else True)
        
        # Otherwise it's a TH1 hopefully
        else:
            titleSize = 0.09
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(kBlack,alpha)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(kBlack,alpha)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)
            
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)

            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.SetMaximum(1.13*hist.GetMaximum())
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.1)
                hist.Draw(datastyle)
                CMS_lumi.CMS_lumi(thisPad, year, 11)
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.35, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.35))

                else:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                if len(signals) == 0:
                    nsignals = 0
                elif addSignals:
                    nsignals = 1
                else:
                    nsignals = len(signals[0])
                legend_topY = 0.73-0.03*(min(len(bkglist[0]),6)+nsignals+1)
                # legend_bottomY = 0.2+0.02*(len(bkglist[0])+nsignals+1)

                legends.append(TLegend(0.65,legend_topY,0.90,0.88))
                legend_duplicates = []
                if not dataOff: legends[hist_index].AddEntry(hist,dataName,datastyle)

                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                if totalBkg == None:
                    tot_hist = hist.Clone(hist.GetName()+'_tot')
                    tot_hist.Reset()
                else:
                    tot_hist = totalBkg[hist_index]

                tot_hist.SetTitle(hist.GetName()+'_tot')
                tot_hist.SetMarkerStyle(0)
                tot_hists.append(tot_hist)
                tot_hists_err.append(tot_hist.Clone())
                tot_hists[hist_index].SetLineColor(kBlack)
                tot_hists_err[hist_index].SetLineColor(kBlack)
                tot_hists_err[hist_index].SetLineWidth(0)
                tot_hists_err[hist_index].SetFillColor(kBlack)
                tot_hists_err[hist_index].SetFillStyle(3354)
                legends[hist_index].AddEntry(tot_hists_err[hist_index],'Total bkg unc.','F')

                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.17)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.1)

                subs[hist_index].SetLeftMargin(0.17)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0)
                subs[hist_index].SetBottomMargin(0.35)
                mains[hist_index].Draw()
                subs[hist_index].Draw()

                # If logy, put QCD on top
                # if logy: bkglist[0], bkglist[-1] = bkglist[-1], bkglist[0]

                # Order based on colors
                if color_idx_order == None:
                    color_idx_order = ColorCodeSortedIndices(colors)
                    colors = [colors[i] for i in color_idx_order]
                
                bkglist[hist_index] = [bkglist[hist_index][i] for i in color_idx_order]
                if bkgNames != [] and isinstance(bkgNames[0],list):
                    bkgNames[hist_index] = [bkgNames[hist_index][i] for i in color_idx_order]

                # Build the stack
                legend_info = collections.OrderedDict()
                for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                    # bkg.Sumw2()
                    if totalBkg == None: tot_hists[hist_index].Add(bkg)
                    
                    if logy:
                        bkg.SetMinimum(1e-3)

                    if 'qcd' in bkg.GetName():
                        bkg.SetFillColor(kYellow)
                        bkg.SetLineColor(kYellow)
                    else:
                        if colors[bkg_index] != None:
                            bkg.SetFillColor(colors[bkg_index])
                            bkg.SetLineColor(colors[bkg_index] if colors[bkg_index]!=0 else kBlack)
                        else:
                            bkg.SetFillColor(default_colors[bkg_index])
                            bkg.SetLineColor(default_colors[bkg_index] if colors[bkg_index]!=0 else kBlack)

                    stacks[hist_index].Add(bkg)
                    if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                    elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                    else: this_bkg_name = bkgNames[hist_index][bkg_index]
                    
                    legend_info[this_bkg_name] = bkg

                # Deal with legend which needs ordering reversed from stack build
                for bname in reversed(legend_info.keys()):
                    if bname not in legend_duplicates:
                        legends[hist_index].AddEntry(legend_info[bname],bname,'f')
                        legend_duplicates.append(bname)
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()

                # Set y max of all hists to be the same to accommodate the tallest
                histList = [stacks[hist_index],tot_hists[hist_index],hist]

                yMax = histList[0].GetMaximum()
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                for h in histList:
                    h.SetMaximum(yMax*(2.5-legend_topY+0.03))
                    if logy == True:
                        h.SetMaximum(yMax*10**(2.5-legend_topY+0.1))

                
                mLS = 0.08
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.15,"xy")
                hist.GetYaxis().SetTitleOffset(1.04)
                hist.GetYaxis().SetLabelSize(0.07)
                hist.GetYaxis().SetTitleSize(titleSize)
                hist.GetXaxis().SetLabelSize(0.07)
                hist.GetXaxis().SetTitleSize(titleSize)
                hist.GetXaxis().SetLabelOffset(0.05)
                if logy == True:
                    hist.SetMinimum(1e-3)
                
                hist.GetYaxis().SetNdivisions(508)

                hist.Draw(datastyle+' X0')
                #gStyle.SetErrorX(0.5)

                if logy == True:stacks[hist_index].SetMinimum(1e-3) 
                
                stacks[hist_index].Draw('same hist') # need to draw twice because the axis doesn't exist for modification until drawing
                try:
                    stacks[hist_index].GetYaxis().SetNdivisions(508)
                except:
                    stacks[hist_index].GetYaxis().SetNdivisions(8,5,0)
                stacks[hist_index].Draw('same hist')

                # Do the signals
                sigs_to_plot = []
                if len(signals) > 0: 
                    # Can add together for total signal
                    if addSignals:
                        totsig = signals[hist_index][0].Clone()
                        for isig in range(1,len(signals[hist_index])):
                            totsig.Add(signals[hist_index][isig])
                        sigs_to_plot = [totsig]
                    # or treat separately
                    else:
                        for sig in signals[hist_index]:
                            sigs_to_plot.append(sig)

                    # Plot either way
                    tot_sigs.append(sigs_to_plot)
                    for isig,sig in enumerate(sigs_to_plot):
                        sig.SetLineColor(kBlack)
                        sig.SetLineWidth(2)
                        if logy == True:
                            sig.SetMinimum(1e-3)
                        # if signalNames == []: this_sig_name = sig.GetName().split('_')[0]
                        if type(signalNames) == str: this_sig_name = signalNames
                        else: this_sig_name = signalNames[isig]

                        legends[hist_index].AddEntry(sig,this_sig_name,'L')
                        sig.Draw('hist same')

                # Draw total hist and error
                if logy: 
                    tot_hists[hist_index].SetMinimum(1e-3)
                    tot_hists_err[hist_index].SetMinimum(1e-3)
                tot_hists[hist_index].Draw('hist same')
                tot_hists_err[hist_index].Draw('e2 same')

                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()

                if not dataOff:
                    #gStyle.SetErrorX(0)
                    hist.Draw(datastyle+' X0 same')
                    #gStyle.SetErrorX(0.5)

                gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
                pulls.append(Make_Pull_plot(hist,tot_hists[hist_index]))
                pulls[hist_index].SetFillColor(kBlue)
                pulls[hist_index].SetTitle(";"+hist.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
                pulls[hist_index].SetStats(0)

                LS = .16

                pulls[hist_index].GetYaxis().SetRangeUser(-2.9,2.9)
                pulls[hist_index].GetYaxis().SetTitleOffset(0.4)
                # pulls[hist_index].GetXaxis().SetTitleOffset(0.9)
                             
                pulls[hist_index].GetYaxis().SetLabelSize(0.13)
                pulls[hist_index].GetYaxis().SetTitleSize(0.12)
                pulls[hist_index].GetYaxis().SetNdivisions(306)
                pulls[hist_index].GetXaxis().SetLabelSize(0.13)
                pulls[hist_index].GetXaxis().SetTitleSize(0.15)

                pulls[hist_index].GetXaxis().SetTitle(xtitle)
                pulls[hist_index].GetYaxis().SetTitle("(Data-Bkg)/#sigma")
                pulls[hist_index].Draw('hist')

                if logy == True:
                    mains[hist_index].SetLogy()

                CMS_lumi.extraText = extraText#'Preliminary'
                CMS_lumi.cmsTextSize = 0.9
                CMS_lumi.cmsTextOffset = 2
                CMS_lumi.lumiTextSize = 0.9
                
                CMS_lumi.CMS_lumi(mains[hist_index], year, 11)
                mains[hist_index].cd()
                lumiE = TLatex()
                lumiE.SetNDC()
                lumiE.SetTextAngle(0)
                lumiE.SetTextColor(kBlack)
                lumiE.SetTextFont(42)
                lumiE.SetTextAlign(31) 
                lumiE.SetTextSize(0.7*0.1)
                lumiE.DrawLatex(1-0.05,1-0.1+0.2*0.1,"137 fb^{-1} (13 TeV)")
                
                if isinstance(subtitles,list) and len(subtitles) > 0:
                    subtitle = TLatex()
                    subtitle.SetNDC()
                    subtitle.SetTextAngle(0)
                    subtitle.SetTextColor(kBlack)
                    subtitle.SetTextFont(42)
                    subtitle.SetTextAlign(12) 
                    subtitle.SetTextSize(0.06)
                    # print (subtitles[hist_index])
                    subtitle_string = '%s < %s < %s %s'%(subtitles[hist_index].split('-')[0], sliceVar.split(' ')[0], subtitles[hist_index].split('-')[1], 'GeV')
                    subtitle.DrawLatex(0.208,0.74,subtitle_string)
    # CMS_lumi.CMS_lumi(myCan, year, 11)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'.pdf','pdf')
        myCan.Print(tag+'/'+name+'.png','png')


def reducedCorrMatrixHist(fit_result,varsOfInterest=[]):
    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if 'Fail_' not in finalPars.at(cm_index).GetName():
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = open('correlation_matrix.txt','w')

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]
        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt.write('%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y]))
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out

def ColorCodeSortedIndices(colors):
    possible_colors = []
    for c in colors:
        if c not in possible_colors:
            possible_colors.append(c)

    new_color_order = []
    for c in possible_colors:
        for idx in [idx for idx,color in enumerate(colors) if color == c]:
            if idx not in new_color_order:
                new_color_order.append(idx)
    return new_color_order

def FindCommonString(string_list):
    to_match = ''   # initialize the string we're looking for/building
    for s in string_list[0]:    # for each character in the first string
        passed = True
        for istring in range(1,len(string_list)):   # compare to_match+s against strings in string_list
            string = string_list[istring]
            if to_match not in string:                  # if in the string, add more
                passed = False
            
        if passed == True:
            to_match+=s

    if to_match[-2] == '_':
        return to_match[:-2] 
    else:
        return to_match[:-1]                # if not, return to_match minus final character

    return to_match[:-2]
        
def Make_Pull_plot( DATA,BKG):
    BKGUP, BKGDOWN = Make_up_down(BKG)
    pull = DATA.Clone(DATA.GetName()+"_pull")
    pull.Add(BKG,-1)
    sigma = 0.0
    FScont = 0.0
    BKGcont = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        FScont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        if FScont>=BKGcont:
            FSerr = DATA.GetBinErrorLow(ibin)
            BKGerr = abs(BKGUP.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FScont<BKGcont:
            FSerr = DATA.GetBinErrorUp(ibin)
            BKGerr = abs(BKGDOWN.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FSerr != None:
            sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
        else:
            sigma = sqrt(BKGerr*BKGerr)
        if FScont == 0.0:
            pull.SetBinContent(ibin, 0.0 )  
        else:
            if sigma != 0 :
                pullcont = (pull.GetBinContent(ibin))/sigma
                pull.SetBinContent(ibin, pullcont)
            else :
                pull.SetBinContent(ibin, 0.0 )
    return pull

def Make_up_down(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down

def checkFitForm(xFitForm,yFitForm):
    if '@0' in xFitForm:
        print 'ERROR: @0 in XFORM. This is reserved for the variable x. Please either replace it with x or start your parameter naming with @1. Quitting...'
        quit()
    if '@0' in yFitForm:
        print 'ERROR: @0 in YFORM. This is reserved for the variable y. Please either replace it with y or start your parameter naming with @1. Quitting...'
        quit()
    if 'x' in yFitForm and 'exp' not in yFitForm:
        print 'ERROR: x in YFORM. Did you mean to write "y"? Quitting...'
        quit()
    if 'y' in xFitForm:
        print 'ERROR: y in XFORM. Did you mean to write "x"? Quitting...'
        quit()

def RFVform2TF1(RFVform,shift=0): # shift tells function how much to shift the indices of the coefficients by
    TF1form = ''
    lookingForDigit = False
    for index,char in enumerate(RFVform):
        # print str(index) + ' ' + char
        if char == '@':
            TF1form+='['
            lookingForDigit = True
        elif lookingForDigit:
            if char.isdigit():
                if RFVform[index+1].isdigit() == False:     # if this is the last digit
                    TF1form+=str(int(char)+shift)
            else:
                TF1form+=']'+char
                lookingForDigit = False
        else:
            TF1form+=char

    return TF1form

# Taken from Kevin's limit_plot_shape.py
def make_smooth_graph(h2,h3):
    h2 = TGraph(h2)
    h3 = TGraph(h3)
    npoints = h3.GetN()
    h3.Set(2*npoints+2)
    for b in range(npoints+2):
        x1, y1 = (ROOT.Double(), ROOT.Double())
        if b == 0:
            h3.GetPoint(npoints-1, x1, y1)
        elif b == 1:
            h2.GetPoint(npoints-b, x1, y1)
        else:
            h2.GetPoint(npoints-b+1, x1, y1)
        h3.SetPoint(npoints+b, x1, y1)
    return h3

def Condor(runscript):
    dir_base = os.environ['CMSSW_BASE']+'/src/2DAlphabet/'

    commands = []

    # tar_cmd = "tar czvf tarball.tgz --directory="+dir_base+" "
    # for t in tarFilesList: tar_cmd+=t+' '
    # commands.append(tar_cmd)

    # Make JDL from template
    timestr = time.strftime("%Y%m%d-%H%M%S")
    out_jdl = 'temp_'+timestr+'_jdl'
    commands.append("sed 's$TEMPSCRIPT$"+runscript+"$g' "+dir_base+"condor/templates/jdl_template_noargs > "+out_jdl)
    commands.append('mkdir notneeded')
    commands.append("condor_submit "+out_jdl)
    commands.append("mv "+out_jdl+" "+dir_base+"condor/jdls/")

    # commands.append("condor_q lcorcodi")

    for s in commands:
        print s
        subprocess.call([s],shell=True)

def StatsForCondor(run_name,toyDict,tarFilesList,commands,files_to_grab=[]):
    if len(files_to_grab) == 0:
        files_to_grab = ['higgsCombine'+run_name+'*.*.mH120.*.root','fitDiagnostics'+run_name+'*.root']
    ntoys = toyDict['ntoys']
    toys_per_job = toyDict['toys_per_job']
    jobs = toyDict['toyjobs']
    seed = toyDict['seed']

    executeCmd('mkdir condor_'+run_name)
    executeCmd('rm condor_'+run_name+'/*')

    with cd('condor_'+run_name):
        # Ship along a post-processing script to grab the output from condor
        tar_cmd = "tar czvf tarball.tgz --directory="+os.environ['CMSSW_BASE']+'/src/2DAlphabet/'+" "
        for t in tarFilesList: tar_cmd+=t+' '
        executeCmd(tar_cmd)

        out_tar_cmd = 'tar -czvf '+run_name+'.tgz '
        for f in files_to_grab: out_tar_cmd+=f+' '

        if jobs > 1:
            for i in range(jobs):
                executeCmd('cp '+os.environ['CMSSW_BASE']+'/src/2DAlphabet/condor/templates/run_combine.sh run_combine_%s.sh'%(i))
                this_run_file = open('run_combine_%s.sh'%(i),'a+')

                this_seed = random.randint(100000,999999)
                this_run_name = run_name+'_%s'%(this_seed)

                these_commands = []
                for cmd in commands:
                    this_command = cmd.replace('-t '+str(ntoys),'-t '+str(toys_per_job))
                    this_command = this_command.replace('-n '+run_name,'-n '+this_run_name)
                    this_command = this_command.replace('-s '+str(seed),'-s '+str(this_seed))
                    # this_command = this_command.replace('-d '+workspace_name,'-d '+projDir+'/'+workspace_name)
                    this_command = this_command.replace('.'+str(seed)+'.','.'+str(this_seed)+'.')
                    this_command = this_command.replace('higgsCombine'+run_name,'higgsCombine'+this_run_name)
                    these_commands.append(this_command)

                this_run_file.write('\n')
                for c in these_commands:
                    this_run_file.write(c+'\n')
                this_run_file.write(out_tar_cmd.replace(run_name+'.tgz',run_name+'_'+str(this_seed)+'.tgz')+'\n')
                this_run_file.write('cp '+this_run_name+'.tgz $CMSSW_BASE/../')
                this_run_file.close()
                Condor('run_combine_%s.sh'%(i))
                
        else:
            executeCmd('cp '+os.environ['CMSSW_BASE']+'/src/2DAlphabet/condor/templates/run_combine.sh run_combine.sh')
            this_run_file = open('run_combine.sh','a+')
            this_run_file.write('\n')
            for c in commands:
                this_run_file.write(c+'\n')
            this_run_file.write(out_tar_cmd+'\n')
            this_run_file.write('cp '+run_name+'.tgz $CMSSW_BASE/../')
            this_run_file.close()
            Condor('run_combine.sh')

def makeShellFinisher(identifier,files_to_grab):
    tar_cmd = 'tar -czvf '+identifier+'.tgz '
    for f in files_to_grab: tar_cmd+=f+' '
    shell_finisher = open('shell_finisher.sh','w')
    shell_finisher.write('#!/bin/sh\n')
    shell_finisher.write('echo "'+tar_cmd+'" \n')
    shell_finisher.write(tar_cmd+'\n')
    shell_finisher.write('cp '+identifier+'.tgz $CMSSW_BASE/../')
    shell_finisher.close()

# Built to wait for condor jobs to finish and then check that they didn't fail
# The script that calls this function will quit if there are any job failures
# listOfJobs input should be whatever comes before '.listOfJobs' for the set of jobs you submitted
def WaitForJobs( listOfJobs ):
    # Runs grep to count the number of jobs - output will have non-digit characters b/c of wc
    preNumberOfJobs = subprocess.check_output('grep "python" '+listOfJobs+' | wc -l', shell=True)
    commentedNumberOfJobs = subprocess.check_output('grep "# python" '+listOfJobs+'.listOfJobs | wc -l', shell=True)

    # Get rid of non-digits and convert to an int
    preNumberOfJobs = int(filter(lambda x: x.isdigit(), preNumberOfJobs))
    commentedNumberOfJobs = int(filter(lambda x: x.isdigit(), commentedNumberOfJobs))
    numberOfJobs = preNumberOfJobs - commentedNumberOfJobs

    finishedJobs = 0
    # Rudementary progress bar
    while finishedJobs < numberOfJobs:
        # Count how many output files there are to see how many jobs finished
        # the `2> null.txt` writes the stderr to null.txt instead of printing it which means
        # you don't have to look at `ls: output_*.log: No such file or directory`
        finishedJobs = subprocess.check_output('ls output_*.log 2> null.txt | wc -l', shell=True)
        finishedJobs = int(filter(lambda x: x.isdigit(), finishedJobs))
        sys.stdout.write('\rProcessing ' + str(listOfJobs) + ' - ')
        # Print the count out as a 'progress bar' that refreshes (via \r)
        sys.stdout.write("%i / %i of jobs finished..." % (finishedJobs,numberOfJobs))
        # Clear the buffer
        sys.stdout.flush()
        # Sleep for one second
        time.sleep(1)


    print 'Jobs completed. Checking for errors...'
    numberOfTracebacks = subprocess.check_output('grep -i "Traceback" output*.log | wc -l', shell=True)
    numberOfSyntax = subprocess.check_output('grep -i "Syntax" output*.log | wc -l', shell=True)

    numberOfTracebacks = int(filter(lambda x: x.isdigit(), numberOfTracebacks))
    numberOfSyntax = int(filter(lambda x: x.isdigit(), numberOfSyntax))

    # Check there are no syntax or traceback errors
    # Future idea - check output file sizes
    if numberOfTracebacks > 0:
        print str(numberOfTracebacks) + ' job(s) failed with traceback error'
        quit()
    elif numberOfSyntax > 0:
        print str(numberOfSyntax) + ' job(s) failed with syntax error'
        quit()
    else:
        print 'No errors!'

def Inter(g1,g2):
    xaxisrange = g1.GetXaxis().GetXmax()-g1.GetXaxis().GetXmin()
    xaxismin = g1.GetXaxis().GetXmin()
    xinters = []
    yinters = []
    for x in range(0,10000):
        xpoint = xaxismin + (float(x)/1000.0)*xaxisrange
        xpoint1 = xaxismin + (float(x+1)/1000.0)*xaxisrange
        Pr1 = g1.Eval(xpoint)
        Pr2 = g2.Eval(xpoint)
        Po1 = g1.Eval(xpoint1)
        Po2 = g2.Eval(xpoint1)
        if (Pr1-Pr2)*(Po1-Po2)<0:
            xinters.append(0.5*(xpoint+xpoint1))
            yinters.append(0.5*(Po1+Po2))
        
    if len(xinters) == 0:
        xinters = [-1]
        yinters = [-1]
    return xinters[0],yinters[0]

# Right after I wrote this I realized it's obsolete... It's cool parentheses parsing though so I'm keeping it
def separateXandYfromFitString(fitForm):
    # Look for all opening and closing parentheses
    openIndexes = []
    closedIndexes = []
    for index,char in enumerate(fitForm):
        if char == '(': # start looking for ")" after this "("
            openIndexes.append(index)
        if char == ')':
            closedIndexes.append(index)


    # Now pair them by looking at the first in closedIndexes and finding the closest opening one (to the left)
    # Remove the pair from the index lists and repeat
    innerContent = []
    for iclose in closedIndexes:
        diff = len(fitForm)     # max length to start because we want to minimize
        for iopen in openIndexes:
            if iclose > iopen:
                this_diff = iclose - iopen
                if this_diff < diff:
                    diff = this_diff
                    candidateOpen = iopen
        openIndexes.remove(candidateOpen)
        innerContent.append(fitForm[iclose-diff+1:iclose])


    outerContent = []
    for c in innerContent:
        keep_c = True
        for d in innerContent:
            if '('+c+')' in d and c != d:
                keep_c = False
                break
        if keep_c:
            outerContent.append(c)

    if len(outerContent) != 2:
        print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
        quit()
    else:
        for c in outerContent:
            if 'x' in c and 'y' not in c:
                xPart = c
            elif 'x' not in c and 'y' in c:
                yPart = c
            else:
                print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
                quit()

    return xPart,yPart

@contextmanager
def cd(newdir):
    print 'cd '+newdir
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
