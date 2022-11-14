from time import time
from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
from TwoDAlphabet.helpers import make_env_tarball, cd, execute_cmd
from TwoDAlphabet.ftest import FstatCalc
import os
import numpy as np

# Helper function to get region names
def _get_other_region_names(pass_reg_name):
    '''
    If passing e.g. "fail", will return ("fail", "pass")
    In HSCP analysis, we're just considering F/P regions
    '''
    return pass_reg_name, pass_reg_name.replace('fail','pass')

# Helper function to generate constraints for parametric Transfer Functions
# Change values as you see fit
def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":-30,"MAX":30}
        else:
            out[i] = {"MIN":-50,"MAX":50}
    return out

# Dict to store transfer function forms and constraints
_rpf_options = {
    '0x0': {
        'form': '0.1*(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '0.1*(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '0.1*(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '0.1*(@0+@1*x)*(1+@2*y)',
        'constraints': _generate_constraints(3)
    },
    '2x0': {
        'form': '0.1*(@0+@1*x+@2*x**2)*(@3)',
        'constraints': _generate_constraints(4)
    },
    '2x1': {
        'form': '0.1*(@0+@1*x+@2*x**2)*(1+@3*y)',
        'constraints': _generate_constraints(4)
    },
    '2x2': {
        'form': '0.1*(@0+@1*x+@2*x**2)*(1+@3*y+@4*y**2)',
        'constraints': _generate_constraints(4)
    }
}

# Helper function for selecting the signal from the ledger
def _select_signal(row, args):
    signame = args[0]
    poly_order = args[1]
    if row.process_type == 'SIGNAL':
        if signame in row.process:
            return True
        else:
            return False
    elif 'Background' in row.process:
        if row.process == 'Background_'+poly_order:
            return True
        elif row.process == 'Background':
            return True
        else:
            return False
    else:
        return True


# Make the workspace
def make_workspace():
    # Create the workspace directory, using info from the specified JSON file
    twoD = TwoDAlphabet('HSCP_fits', 'test.json', loadPrevious=False)

    # 2DAlphabet wasn't intended for an analysis like this, so the default function 
    # for Looping over all regions and for a given region's data histogram, subtracting
    # the list of background histograms, and returning a data-bkgList is called initQCDHists.
    # This is b/c QCD multijet is the main background we usually estimate via 2DAlphabet
    bkg_hists = twoD.InitQCDHists()
    print('bkg_hists = {}'.format(bkg_hists))

    # Now, we loop over "pass" and "fail" regions and get the binnings
    for f, p in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'fail' in r]:
	print(f, p)
	# get the binning for the fail region
	binning_f, _ = twoD.GetBinningFor(f)
	# you can change the name as you see fit 
	fail_name = 'Background_'+f
	# this is the actual binned distribution of the fail
	bkg_f = BinnedDistribution(fail_name, bkg_hists[f], binning_f, constant=False)
	# now we add it to the 2DAlphabet ledger
	twoD.AddAlphaObj('Background',f, bkg_f)

        # now construct all of the possible transfer functions, to be chosen and used later
        for opt_name, opt in _rpf_options.items():
           bkg_rpf = ParametricFunction(
               fail_name.replace('fail','rpf')+'_'+opt_name, # this is our pass/fail ratio
               binning_f,                                    # we use the binning from fail
               opt['form'],                                  # was _rpf_options['0x0']['form'],
               constraints=opt['constraints']                # was _rpf_options['0x0']['constraints'] 
           )


           # now define the bkg in pass as the bkg in fail multiplied by the transfer function (bkg_rpf)
           bkg_p = bkg_f.Multiply(fail_name.replace('fail','pass')+'_'+opt_name, bkg_rpf)

           # then add this to the 2DAlphabet ledger
           twoD.AddAlphaObj('Background_'+opt_name,p,bkg_p,title='Background')

    # and save it out
    twoD.Save()


# function for perfomring the fit
def perform_fit(signal, tf, extra=''):
    '''
        signal [str] = 'Type-Mass'
        tf [str] = 0x0, 0x1, 1x0, 1x1, 1x2, 2x2
	extra (str) = any extra flags to pass to Combine when running the ML fit
    '''
    # this is the name of the directory created in the workspace function
    working_area = 'HSCP_fits'
    # we reuse the workspace from the last step.
    # The runConfig.json is copied from the origin JSON config file, 
    # and we must specify that we want to load the previous workspace
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)

    # Now we create a ledger and make a new area for it with a Combine card
    # this select() method uses lambda functions. Will explain later
    subset = twoD.ledger.select(_select_signal, 'HSCP-{}'.format(signal), tf)   
    twoD.MakeCard(subset, 'HSCP-{}-{}_area'.format(signal, tf))

    # perform fit
    print("perform fit")
    #twoD.MLfit('HSCP-{}-{}_area'.format(signal, tf), rMin=0, rMax=30, verbosity=1, extra=extra)
    twoD.MLfit('HSCP-{}-{}_area'.format(signal, tf), rMin=0, rMax=10, verbosity=1, extra=extra)

def plot_fit(signal, tf):
    working_area = 'HSCP_fits'
    print("DoingTwoDAlphabet")
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    print("Doing twoD.ledger.select")
    subset = twoD.ledger.select(_select_signal, 'HSCP-{}'.format(signal), tf) 
    print("Doing twoD.StdPlots")
    twoD.StdPlots('HSCP-{}-{}_area'.format(signal, tf), subset)
    twoD.StdPlots('HSCP-{}-{}_area'.format(signal, tf), subset, True)

def GOF(signal,tf,condor=True, extra=''):
    # replace the blindedFit option in the config file with COMMENT to effectively "unblind" the GoF
    findReplace = {"blindedFit": "COMMENT"}
    working_area = 'HSCP_fits'
    signame = 'HSCP-'+signal
    #twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area),findreplace=findReplace, loadPrevious=True)
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    if not os.path.exists(twoD.tag+'/'+signame+'-{}_area/card.txt'.format(tf)):
        print('{}/{}-area/card.txt does not exist, making card'.format(twoD.tag,signame))
        subset = twoD.ledger.select(_select_signal, signame, tf)
        twoD.MakeCard(subset, signame+'_area')
    if condor == False:
        twoD.GoodnessOfFit(
            signame+'-{}_area'.format(tf), ntoys=500, freezeSignal=0,
            condor=False
        )
    else:
        twoD.GoodnessOfFit(
            signame+'-{}_area'.format(tf), ntoys=500, freezeSignal=0,
            condor=True, njobs=10
        )

def plot_GOF(condor=True):
    working_area = 'HSCP_fits'
    plot.plot_gof('{}'.format(working_area), 'HSCP-{}-{}_area'.format(signal, tf), condor=condor)

def load_RPF(twoD):
    '''	
	loads the rpf parameter values for use in toy generation
    '''
    params_to_set = twoD.GetParamsOnMatch('rpf.*', 'Signal', 'b')
    return {k:v['val'] for k,v in params_to_set.items()}

def SignalInjection(r, condor=True):
    working_area = 'HSCP_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    #params = load_RPF(twoD)
    twoD.SignalInjection(
	'HSCP-{}-{}_area'.format(signal, tf),
	injectAmount = r,	# injected signal xsec (r=0 : bias test)
	ntoys=500,		# will take forever if not on condor
	blindData = True,	# make sure you're blinding if working with data
	#setParams = params,     # give the toys the same RPF params
	verbosity = 0,		# you can change this if you need
	njobs=10,
	condor = condor
    )

def plot_SignalInjection(r, condor=False):
    working_area = 'HSCP_fits'
    plot.plot_signalInjection(working_area, 'HSCP-{}-{}_area'.format(signal, tf), injectedAmount=r, condor=condor)

def Impacts():
    working_area = 'HSCP_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    twoD.Impacts('HSCP-{}-{}_area'.format(signal, tf), cardOrW='card.txt', extra='-t 1')

def run_limits():
    working_area = 'HSCP_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    twoD.Limit(
	subtag='HSCP-{}-{}_area'.format(signal, tf),
	blindData=False,	# BE SURE TO CHANGE THIS IF YOU NEED TO BLIND YOUR DATA 
	verbosity=1,
	condor=False
    )

def _gof_for_FTest(twoD, subtag, card_or_w='card.txt'):

    run_dir = twoD.tag+'/'+subtag
    
    with cd(run_dir):
        gof_data_cmd = [
            'combine -M GoodnessOfFit',
            '-d '+card_or_w,
            '--algo=saturated',
            '-n _gof_data'
        ]

        gof_data_cmd = ' '.join(gof_data_cmd)
        execute_cmd(gof_data_cmd)

def test_FTest(poly1, poly2, signal=''):
    '''
    Perform an F-test using existing working areas
    '''
    working_area = 'HSCP_fits'
    
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    binning = twoD.binnings['default']
    nBins = (len(binning.xbinList)-1)*(len(binning.ybinList)-1)
    
    # Get number of RPF params and run GoF for poly1
    params1 = twoD.ledger.select(_select_signal, 'HSCP-{}'.format(signal), poly1).alphaParams
    rpfSet1 = params1[params1["name"].str.contains("rpf")]
    print("rpfSet1: " + str(rpfSet1))
    nRpfs1  = len(rpfSet1.index)
    print(" >>>>>> Num RPF parameters for poly1: " + str(nRpfs1))
    _gof_for_FTest(twoD, 'HSCP-gluino-1800-{}_area'.format(poly1), card_or_w='card.txt')
    gofFile1 = working_area+'/HSCP-gluino-1800-125-{}_area/higgsCombine_gof_data.GoodnessOfFit.mH120.root'.format(poly1)

    # Get number of RPF params and run GoF for poly2
    params2 = twoD.ledger.select(_select_signal, 'HSCP-{}'.format(signal), poly2).alphaParams
    rpfSet2 = params2[params2["name"].str.contains("rpf")]
    nRpfs2  = len(rpfSet2.index)
    print(" >>>>>> Num RPF parameters for poly2: " + str(nRpfs2))
    _gof_for_FTest(twoD, 'HSCP-gluino-1800-{}_area'.format(poly2), card_or_w='card.txt')
    gofFile2 = working_area+'/HSCP-gluino-1800-{}_area/higgsCombine_gof_data.GoodnessOfFit.mH120.root'.format(poly2)


    base_fstat = FstatCalc(gofFile1,gofFile2,nRpfs1,nRpfs2,nBins)
    print(base_fstat)

    def plot_FTest(base_fstat,nRpfs1,nRpfs2,nBins):
        from ROOT import TF1, TH1F, TLegend, TPaveText, TLatex, TArrow, TCanvas, kBlue, gStyle
        gStyle.SetOptStat(0000)

        if len(base_fstat) == 0: base_fstat = [0.0]

        ftest_p1    = min(nRpfs1,nRpfs2)
        ftest_p2    = max(nRpfs1,nRpfs2)
        ftest_nbins = nBins
        fdist       = TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(10,1.3*base_fstat[0]))
        fdist.SetParameter(0,1)
        fdist.SetParameter(1,ftest_p2-ftest_p1)
        fdist.SetParameter(2,ftest_nbins-ftest_p2)

        pval = fdist.Integral(0.0,base_fstat[0])
        print 'P-value: %s'%pval

        c = TCanvas('c','c',800,600)    
        c.SetLeftMargin(0.12) 
        c.SetBottomMargin(0.12)
        c.SetRightMargin(0.1)
        c.SetTopMargin(0.1)
        ftestHist_nbins = 30
        ftestHist = TH1F("Fhist","",ftestHist_nbins,0,max(10,1.3*base_fstat[0]))
        ftestHist.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
        ftestHist.GetXaxis().SetTitleSize(0.025)
        ftestHist.GetXaxis().SetTitleOffset(2)
        ftestHist.GetYaxis().SetTitleOffset(0.85)
        
        ftestHist.Draw("pez")
        ftestobs  = TArrow(base_fstat[0],0.25,base_fstat[0],0)
        ftestobs.SetLineColor(kBlue+1)
        ftestobs.SetLineWidth(2)
        fdist.Draw('same')

        ftestobs.Draw()
        tLeg = TLegend(0.6,0.73,0.89,0.89)
        tLeg.SetLineWidth(0)
        tLeg.SetFillStyle(0)
        tLeg.SetTextFont(42)
        tLeg.SetTextSize(0.03)
        tLeg.AddEntry(ftestobs,"observed = %.3f"%base_fstat[0],"l")
        tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
        tLeg.Draw("same")

        model_info = TPaveText(0.2,0.6,0.4,0.8,"brNDC")
        model_info.AddText('p1 = '+poly1)
        model_info.AddText('p2 = '+poly2)
        model_info.AddText("p-value = %.2f"%(1-pval))
        model_info.Draw('same')
        
        latex = TLatex()
        latex.SetTextAlign(11)
        latex.SetTextSize(0.06)
        latex.SetTextFont(62)
        latex.SetNDC()
        latex.DrawLatex(0.12,0.91,"CMS")
        latex.SetTextSize(0.05)
        latex.SetTextFont(52)
        latex.DrawLatex(0.2,0.91,"Preliminary")
        latex.SetTextFont(42)
        latex.SetTextFont(52)
        latex.SetTextSize(0.045)
        c.SaveAs(working_area+'/ftest_{0}_vs_{1}_notoys.png'.format(poly1,poly2))

    plot_FTest(base_fstat,nRpfs1,nRpfs2,nBins)

if __name__ == "__main__":
    make_workspace()
    perform_fit('gluino-1800','0x0',extra='--robustHesse 1')
    #plot_fit('gluino-1800','0x0' )
    #GOF(condor=False, extra='--text2workspace "--channel-masks" --setParametersForFit mask_pass_SIG=1 --setParametersForEval mask_pass_SIG=1')
    #GOF(condor=False, extra='')
    #SignalInjection(0, condor=False)	# you can make a loop to run a bunch of injected xsecs
    #run_limits()
    #Impacts()

    #test_FTest('0x0','1x0')


    # if you ran GOF/SigInj via condor, you need to wait until they're finished to run plotting:
    #plot_GOF(condor=False)
    #plot_SignalInjection(0, condor=False)
