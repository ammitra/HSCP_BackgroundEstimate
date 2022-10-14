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

    # Now let's make the parametric transfer function to go from Fail -> Pass
    bkg_rpf = ParametricFunction(
	fail_name.replace('fail','rpf'),	# this is our pass/fail ratio
	binning_f,				# we use the binning from fail
	_rpf_options['0x0']['form'],		# let's make it constant in Ias (and ProbQ obviously)
	_rpf_options['0x0']['constraints']	# use the default constraints [0,5]
    )

    # now define the bkg in pass as the bkg in fail multiplied by the transfer function (bkg_rpf)
    bkg_p = bkg_f.Multiply(fail_name.replace('fail','pass'), bkg_rpf)

    # then add this to the 2DAlphabet ledger
    twoD.AddAlphaObj('Background',p,bkg_p,title='Background')

    # and save it out
    twoD.Save()


# function for perfomring the fit
def perform_fit(extra=''):
    '''
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
    subset = twoD.ledger.select(_select_signal, 'Signal', '')   
    twoD.MakeCard(subset, 'Signal_area')

    # perform fit
    twoD.MLfit('Signal_area', rMin=-1, rMax=20, verbosity=1, extra=extra)

def plot_fit():
    working_area = 'HSCP_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, 'Signal', '')
    twoD.StdPlots('Signal_area', subset)

if __name__ == "__main__":
    make_workspace()
    perform_fit(extra='--robustHesse 1 --freezeParameters var{Background_fail.*}')
    plot_fit()
