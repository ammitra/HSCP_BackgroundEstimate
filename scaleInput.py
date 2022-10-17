import ROOT

def addEvents(hist):
    for i in range(1,21):
	if hist.GetBinContent(i, 1) == 0:
	    hist.SetBinContent(i, 1, 20)

def Scale(amount):
    for proc in ['Data']:
	fTemp = ROOT.TFile.Open('HSPC_{}.root'.format(proc),'READ')
	print(fTemp)

	hFailTemp = fTemp.Get('hfail')
	print(hFailTemp)
	hFail = hFailTemp.Clone('hfail')
	print(hFail)
	hFail.Scale(amount)
	hFail.SetDirectory(0)
	print(hFail)

	hPassTemp = fTemp.Get('hpass')
	print(hPassTemp)
	hPass = hPassTemp.Clone('hpass')
	print(hPass)
	hPass.Scale(amount)
	hPass.SetDirectory(0)
	print(hPass)

	fTemp.Close()	

	fNew = ROOT.TFile.Open('HSPC_{}_scaled.root'.format(proc),'RECREATE')

	addEvents(hPass)
	addEvents(hFail)

	hPass.Write()
	hFail.Write()
	fNew.Close()

Scale(10.0)
