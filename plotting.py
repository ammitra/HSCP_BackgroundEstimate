def plotRPF(postfitShapesFile,odir):
    hPass = get2DPostfitPlot(postfitShapesFile,"Background","pass")
    hFail = get2DPostfitPlot(postfitShapesFile,"Background","fail")
    hPass = hPass.ProjectionX("hPass_for_rpf")
    hFail = hFail.ProjectionX("hFail_for_rpf")

    hRpf  = hPass.Clone("hRPF")
    hRpf.Divide(hFail)

    hRpf, edges = hist2array(hRpf,return_edges=True)
    #print(hRpf, edges)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    hep.histplot(hRpf,edges[0],label = "$R_{P/F}$", histtype="step",edgecolor='blue')

    hep.cms.text("Work in progress",loc=0)
    ax.set_xlim([0.,1.])
    ax.set_ylim([0,0.5])


    plt.xlabel("$I_{as}$",horizontalalignment='right', x=1.0)
    plt.ylabel("$R_{P/F}$",horizontalalignment='right', y=1.0)
    ax.yaxis.set_tick_params(which='minor', left=False)    
    ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"/RPF.pdf")
    plt.savefig(odir+"/RPF.pdf", bbox_inches='tight')
    plt.savefig(odir+"/RPF.png", bbox_inches='tight')
    plt.cla()
    plt.clf()

def getPoissonErrors(hist,binWidthDivision=False):
    hist.SetBinErrorOption(1)

    #This is needed due to some nasty float precision inaccuracy causing some data content to be 0.9999998
    #The kPoisson error does not get calculated correctly in that case for some reason
    tempHist   = hist.Clone("tempHist_forErrs")
    tempHist.Reset()
    tempHist.SetBinErrorOption(1)

    errors_low = []
    errors_hi  = []
    for i in range(1,hist.GetNbinsX()+1):
        tempHist.SetBinContent(i,int(round(hist.GetBinContent(i))))
        #print(int(hist.GetBinContent(i)),tempHist.GetBinErrorLow(i),tempHist.GetBinErrorUp(i))
        if(binWidthDivision):
            errors_low.append(tempHist.GetBinErrorLow(i)/tempHist.GetBinWidth(i))
            errors_hi.append(tempHist.GetBinErrorUp(i)/tempHist.GetBinWidth(i))
        else:
            errors_low.append(tempHist.GetBinErrorLow(i))
            errors_hi.append(tempHist.GetBinErrorUp(i))

    return [errors_low,errors_hi]

def rebinHisto(hModel,hToRebin,name,scale=1.0):
    hRes = hModel.Clone(name)
    hRes.Reset()
    xaxis = hToRebin.GetXaxis()
    yaxis = hToRebin.GetYaxis()
    xaxis_re = hRes.GetXaxis()
    yaxis_re = hRes.GetYaxis()
    for i in range(1,hToRebin.GetNbinsX()+1):
        for j in range(1,hToRebin.GetNbinsY()+1):
            x = xaxis.GetBinCenter(i)
            y = yaxis.GetBinCenter(j)
            i_re = xaxis_re.FindBin(x)
            j_re = yaxis_re.FindBin(y)
            value = hToRebin.GetBinContent(i,j)
            if(value<0.):
                value = 0.
            err = hToRebin.GetBinError(i,j)
            err_re = np.sqrt(hRes.GetBinError(i_re,j_re)*hRes.GetBinError(i_re,j_re)+err*err)
            hRes.Fill(x,y,value)
            hRes.SetBinError(i_re,j_re,err_re)
    hRes.Scale(scale)
    hRes.SetDirectory(0)
    return hRes

def get_binning_x(hLow,hSig,hHigh):
    bins = []
    for i in range(1,hLow.GetNbinsX()+1):
        bins.append(hLow.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hSig.GetNbinsX()+1):
        bins.append(hSig.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hHigh.GetNbinsX()+2):#low edge of overflow is high edge of last bin
        bins.append(hHigh.GetXaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def get_binning_y(hLow,hSig,hHigh):
    #histos should have same binning in Y
    bins = []
    for i in range(1,hLow.GetNbinsY()+2):
        bins.append(hLow.GetYaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def merge_low_sig_high(hLow,hSig,hHigh,hName="temp"):
    n_x_low     = hLow.GetNbinsX()
    n_x_sig     = hSig.GetNbinsX()
    n_x_high    = hHigh.GetNbinsX()
    n_x         = n_x_low + n_x_sig + n_x_high
    n_y         = hLow.GetNbinsY()#assumes Y bins are the same
    bins_x      = get_binning_x(hLow,hSig,hHigh)
    bins_y      = get_binning_y(hLow,hSig,hHigh)
    h_res       = r.TH2F(hName,"",n_x,bins_x,n_y,bins_y)
    for i in range(1,n_x_low+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+0,j,hLow.GetBinContent(i,j))
            h_res.SetBinError(i+0,j,hLow.GetBinError(i,j))

    for i in range(1,n_x_sig+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_low,j,hSig.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_low,j,hSig.GetBinError(i,j))

    for i in range(1,n_x_high+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_sig+n_x_low,j,hHigh.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_sig+n_x_low,j,hHigh.GetBinError(i,j))
    return h_res

def get2DPostfitPlot(file,process,region):
    f       = r.TFile.Open(file)
    hLow    = f.Get("{0}_LOW_postfit/{1}".format(region,process))
    hSig    = f.Get("{0}_SIG_postfit/{1}".format(region,process))
    hHigh   = f.Get("{0}_HIGH_postfit/{1}".format(region,process))
    h2      = merge_low_sig_high(hLow,hSig,hHigh,hName="h2_{0}_{1}".format(process,region))
    h2.SetDirectory(0)
    return h2

def get2DPrefitPlot(file,process,region):
    f       = r.TFile.Open(file)
    hLow    = f.Get("{0}_LOW_prefit/{1}".format(region,process))
    hSig    = f.Get("{0}_SIG_prefit/{1}".format(region,process))
    hHigh   = f.Get("{0}_HIGH_prefit/{1}".format(region,process))
    h2      = merge_low_sig_high(hLow,hSig,hHigh,hName="h2_{0}_{1}".format(process,region))
    h2.SetDirectory(0)
    return h2

def getUncBand(totalHistos):
    yLo = []
    yUp = []
    for i in range(1,totalHistos.GetNbinsX()+1):
        errLo  = totalHistos.GetBinErrorLow(i)
        errUp  = totalHistos.GetBinErrorUp(i)
        mcPred = totalHistos.GetBinContent(i)
        yLo.append(mcPred-errLo)
        yUp.append(mcPred+errUp)
    return np.array(yLo), np.array(yUp)

def calculatePull(hData,dataErrors,hTotBkg,uncBand):
    pulls = []
    for i,dataYield in enumerate(hData):
        mcYield     = hTotBkg[i]
        diff        = dataYield-mcYield
        dataErr     = np.sqrt(dataYield)
        if(dataYield>=mcYield):
            dataErr = dataErrors[0][i]#ErrorLo
            mcErr   = uncBand[1][i]-mcYield#ErrorUp
        else:
            dataErr = dataErrors[1][i]#ErrorUp
            mcErr   = uncBand[0][i]-mcYield#ErrorLo

        sigma =  np.sqrt(dataErr*dataErr+mcErr*mcErr)
        pull        = diff/sigma
        pulls.append(pull)

    return np.array(pulls)

def divideByBinWidth(h):
    hist, edges = hist2array(h,return_edges=True)
    newHist     = []
    for i in range(len(hist)):
        binWidth    = edges[0][i+1]-edges[0][i]
        newHist.append(hist[i]/binWidth)

    return newHist,edges



def plotShapes(hData,hMC,uncBand,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",logFlag=False):
    dataErrors      = getPoissonErrors(hData,binWidthDivision=True)
    hData, edges    = divideByBinWidth(hData)
    centresData     = (edges[0][:-1] + edges[0][1:])/2.#Bin centers
    xerrorsData     = []

    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)

    histosMC        = []
    for h in hMC:
        histosMC.append(divideByBinWidth(h)[0])


    plt.style.use([hep.style.CMS])
    #matplotlib.rcParams.update({'font.size': 30})
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.05})
    axs = axs.flatten()
    plt.sca(axs[0])

    hep.histplot(histosMC[:-1],edges[0],stack=True,ax=axs[0],label = labelsMC, histtype="fill",facecolor=colorsMC,edgecolor='black')
    plt.errorbar(centresData,hData, yerr=dataErrors, xerr=xerrorsData, fmt='o',color="k",label = "Data")

    for i in range(len(uncBand[0])):
        binWidth            = edges[0][i+1]-edges[0][i]
        uncBand[0][i]       = uncBand[0][i]/binWidth
        uncBand[1][i]       = uncBand[1][i]/binWidth

    uncBandLow = np.append(uncBand[0],[0],axis=0)
    uncBandHi  = np.append(uncBand[1],[0],axis=0)#Hack to get the last bin uncertainty to plot, since we're using step="post"

    plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post")

    axs[0].legend()
    plt.ylabel("Events/bin width",horizontalalignment='right', y=1.0)
    axs[1].set_ylabel("Pulls")

    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        if(logFlag):
            yMaximum = max(hData)*10.
            axs[0].set_ylim([1.,yMaximum])
        else:
            yMaximum = max(hData)*1.5+10.
            axs[0].set_ylim([0.,yMaximum])

    if(logFlag):
        axs[0].set_yscale("log")

    if("16APV" in outputFile):
        lumi = 19.5
    elif("16" in outputFile):
        lumi = 16.8
    elif("17" in outputFile):
        lumi = 41.5
    elif("18" in outputFile):
        lumi = 59.8
    else:
        lumi = 138

    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc=1,ncol=2)

    if(projectionText):
        plt.text(0.60, 0.15, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)


    plt.sca(axs[1])#switch to lower pad
    #axs[1].axhline(y=0.0, xmin=0, xmax=1, color="r")
    axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=-1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].axhline(y=-2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].set_ylim([-2.5,2.5])
    plt.xlabel(xlabel,horizontalalignment='right', x=1.0)

    pulls = calculatePull(hData,dataErrors,histosMC[-1],uncBand)
    hep.histplot(pulls,edges[0],ax=axs[1],linewidth=1,histtype="fill",facecolor="grey",edgecolor='black')

    print("Saving ", outputFile)
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))

    plt.clf()
    plt.cla()

def printYields(data_proj,hMC,procs):
    for i,proc in enumerate(procs):
        yieldErr  = ctypes.c_double(1.)#ROOT thing
        procYield = hMC[i].IntegralAndError(1,hMC[i].GetNbinsX(),yieldErr,"")
        print("{0} & {1:.0f} & {2:.0f}".format(proc,procYield,yieldErr.value))
    
    yieldErr  = ctypes.c_double(1.)#ROOT thing
    procYield = data_proj.IntegralAndError(1,data_proj.GetNbinsX(),yieldErr,"")
    print("Data & {0:.0f} & {1:.0f}".format(procYield,yieldErr.value))

def plotPostfit(postfitShapesFile,region,odir):

    labels              = ["Data","Background","Signal"]
    tags                = ["data_obs","Background","Signal"]
    colors              = ["black","khaki","blueviolet"]
    plotSlices          = True



    twoDShapes          = []
    #Exclude signal from considered processes in b-only fit
    if("_b.root" in postfitShapesFile):
        totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalBkg",region)
        labels      = labels[:-1]
        tags        = tags[:-1]
        colors      = colors[:-1]
    else:
        totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalProcs",region)


    #Merge sliced histograms
    for tag in tags:
        twoDShape = get2DPostfitPlot(postfitShapesFile,tag,region)
        twoDShapes.append(twoDShape)


    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        projections.append(proj)
    totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    projections.append(totalProcs_proj)
    uncBand_proj        = getUncBand(totalProcs_proj)

    plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$I_{as}$ [GeV]","{0}/postfit_{1}.png".format(odir,region),logFlag=True)

    tags.append("Total")
    print("Yields in {0}".format(region))
    printYields(projections[0],projections[1:],tags[1:])


if __name__ == '__main__':
    from pathlib import Path
    import matplotlib
    import ROOT as r
    import mplhep as hep
    import numpy as np
    from root_numpy import hist2array
    import matplotlib.pyplot as plt
    import ctypes

    r.gROOT.SetBatch(True)
    matplotlib.use('Agg')
  
    fitDir      = "HSCP_fits/Signal_area/"
    Path(fitDir).mkdir(parents=True, exist_ok=True)

    #b-only plots
    fitFile      = fitDir+"postfitshapes_b.root"
    oDir         = fitDir + "/plots_fit_b/"
    Path(oDir).mkdir(parents=True, exist_ok=True)
    plotPostfit(fitFile,"pass",oDir)
    plotPostfit(fitFile,"fail",oDir)
    plotRPF(fitFile,oDir)

    #s+b plots
    fitFile      = fitDir+"postfitshapes_s.root"
    oDir         = fitDir + "/plots_fit_s/"
    Path(oDir).mkdir(parents=True, exist_ok=True)
    plotPostfit(fitFile,"pass",oDir)    
    plotPostfit(fitFile,"fail",oDir)
    plotRPF(fitFile,oDir)


