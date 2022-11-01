# HSCP Background Estimate

## Getting Started
First, install [2DAlphabet](https://github.com/ammitra/2DAlphabet) from my fork, with the `binningFloatFix` branch that addresses some of the issues with using non-integer bin widths.

Log on to the LPC and ensure that your Github SSH keys are added to the agent.
```
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_XYZ
```

Then, install 2DAlphabet following the modified [instructions](https://github.com/ammitra/2DAlphabet#installation), shown here:
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src
cmsenv
git clone --branch binningFloatFix https://github.com/ammitra/2DAlphabet.git
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
curl -s https://raw.githubusercontent.com/lcorcodilos/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh | bash
scram b clean; scram b -j 4
cmsenv
```

Finally, create the virtual environment into which we install 2DAlphabet:
```
python -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet
python setup.py develop
```

Now, you are ready to use 2DAlphabet. Note that, every time you log on to the LPC you must run these steps to use 2DAlphabet:
```
cd /path/to/CMSSW_10_6_14/src
cmsenv
source twoD-env/bin/activate
```

## Making the template histograms
The HSCP input histograms are created from the templates found in the `crab_Analysis_2018_*.root` files. The script `makeRegions.py` will find the appropriate histogram from the background (`crab_Analysis_2018_AllBackground_CodeV40p4_v1.root`) and signal (`crab_Analysis_2018_HSCPgluino_M-1800_CodeV40p4_v1.root`) files and create new template files for 2DAlphabet. 

The script creates two new ROOT files - one for signal and one for data (we just treat the background as data for this example). Each file contains an `hpass` and `hfail` histogram, where the Pass region is defined as `0 < ProbQ < 0.1` and the Fail as `0.1 < ProbQ < 0.7`. 

The Pass and Fail histograms are 20x1 (`I_as` vs `ProbQ`), where I've moved `I_as` to the X-axis because at the moment 2DAlphabet only allows blinding on the X-axis. 

## 2DAlphabet configuration file
The 2DAlphabet config file is `test.json`. The binning is defined such that `LOW`, `SIG` and `HIGH` each have events in them (otherwise the fit will fail). **2DAlphabet should be modified so that we can pass an empty dict to `HIGH`**

More explanations to come 

## 2DAlphabet python API
To make the workspace and run the fit, run 
```
python test.py
```
explanation to come

## Plotting environment
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0/src
cmsenv
python3 -m virtualenv plotting-env
source plotting-env/bin/activate
pip install numpy root_numpy mplhep matplotlib
```


### Setting up plotting environment after new login
```
cd CMSSW_12_3_0/src
cmsenv
source plotting-env/bin/activate
```
