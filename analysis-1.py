#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

##################################################
import ROOT
import math
#------------------------------------------------#
# CONFIGURATION
#------------------------------------------------#
input_file = ROOT.TFile("montecarlo.root")
num_events = -1 #10000
#------------------------------------------------#
simul = input_file.Get("simul")
truth = input_file.Get("truth")

max_events = simul.GetEntries()

if (num_events < 0) or (num_events > max_events):
    num_events = max_events

print "number of events:", num_events

## the events in the simul tree and the truth tree
## are out of sync. ie. the nth entry in the simul
## tree and the nth entry in the truth tree do not
## necessarily correspond to the same event.
## the line below allows us to find the 
truth.BuildIndex("runNumber","eventNumber")
#------------------------------------------------#
# BOOK HISTOGRAMS
#------------------------------------------------#
h_mu_eta = ROOT.TH1F("h_mu_eta","",100,-5.0,5.0)
h_mu_phi = ROOT.TH1F("h_mu_phi","",100,-3.2,3.2)
pixel_mu = ROOT.TH2F("pixel_mu","",100,-5.0,5.0,100,-3.2,3.2)
pixel_el = ROOT.TH2F("pixel_el","",100,-5.0,5.0,100,-3.2,3.2)
cone_rad_mu = ROOT.TH1F("cone_rad_mu","",1000,0.0,1.0)
cone_rad_el = ROOT.TH1F("cone_rad_el","",1000,0.0,1.0)
mu_ptbins = []
frac_pT_mu  = ROOT.TH2F("frac_pT_mu","",40,20.0,100.,100,0.0,0.2)
frac_pT_el  = ROOT.TH2F("frac_pT_el","",40,20.0,100.,100,0.0,0.2)
#------------------------------------------------#
# ANALYSE DATA
#------------------------------------------------#
for n in range(num_events):
    
    if (n % (num_events/100)) == 0:
        print "done... %3i%%" % (n*100/num_events)

    ## load data from nth entry in the simul tree
    simul.GetEntry(n)
    runNumber   = simul.runNumber
    eventNumber = simul.eventNumber

    ## load data from the entry in the truth tree
    ## that has the same event number and run nu-
     ## mber ie. the entry that corresponds to the
    ## same event
    truth.GetEntryWithIndex(runNumber,eventNumber)

    ## sanity checks
    if runNumber   != truth.runNumber:
        continue
    if eventNumber != truth.eventNumber:
        continue
    #--------------------------------------------#
    ##MUONS
    for i in range(len(simul.mu_pt)):
        
        h_mu_eta.Fill(simul.mu_eta[i])
        h_mu_phi.Fill(simul.mu_phi[i])

        pixel_mu.Fill(simul.mu_eta[i],simul.mu_phi[i])

        for j in range(len(truth.mu_pt)):
            R2_mu=(simul.mu_eta[i]-truth.mu_eta[j])**2+(simul.mu_phi[i]-truth.mu_phi[j])**2
            R_mu=math.sqrt(R2_mu)
            cone_rad_mu.Fill(R_mu)
            if cone_rad_mu>0.1:
                continue
            frac_mu=abs(truth.mu_pt[j]-simul.mu_pt[i])/truth.mu_pt[j]
            mu_pt=truth.mu_pt[j]/1000.
            frac_pT_mu.Fill(mu_pt,frac_mu)

    ##ELECTRONS
    for i in range(len(simul.el_pt)):

        pixel_el.Fill(simul.el_eta[i],simul.el_phi[i])

        for j in range(len(truth.el_pt)):
            R2_el=(simul.el_eta[i]-truth.el_eta[j])**2+(simul.el_phi[i]-truth.el_phi[j])**2
            R_el=math.sqrt(R2_el)
            cone_rad_el.Fill(R_el)
            if cone_rad_el>0.1:
                continue
            frac_el=abs(truth.el_pt[j]-simul.el_pt[i])/truth.el_pt[j]
            el_pt=truth.el_pt[j]/1000.
            frac_pT_el.Fill(el_pt,frac_el)

##MUONS 
nbinsx_mu = frac_pT_mu.GetXaxis().GetNbins()
nbinsy_mu = frac_pT_mu.GetYaxis().GetNbins()
wsums_mu = []

for binx in range(nbinsx_mu):
    wsum=0
    for biny in range(nbinsy_mu):
        wsum +=frac_pT_mu.GetBinContent(binx,biny)
    wsums_mu.append(wsum)
x_pT_mu = []
y_frac_mu = []
error_mu = []
for binx in range(nbinsx_mu):
    mean=0
    mean_squared=0
    if wsums_mu[binx]==0.0:
        continue
    for biny in range(nbinsy_mu):        
        bincenter_y = frac_pT_mu.GetYaxis().GetBinCenter(biny)
        weight = frac_pT_mu.GetBinContent(binx,biny)/wsums_mu[binx]
        mean += weight*bincenter_y
        mean_squared +=weight*bincenter_y**2 
    y_frac_mu.append(mean)
    x_pT_mu.append(frac_pT_mu.GetXaxis().GetBinCenter(binx))
    error_mu.append(math.sqrt(mean_squared-mean**2)/math.sqrt(wsums_mu[binx])) #standard error

mu_graph = ROOT.TGraphErrors(len(x_pT_mu))
mu_graph.SetName("mu_graph")
assert len(x_pT_mu) == len(y_frac_mu)
for p in range(len(x_pT_mu)):
    mu_graph.SetPoint(p, x_pT_mu[p], y_frac_mu[p])
    mu_graph.SetPointError(p, 0, error_mu[p])

##ELECTRONS
nbinsx_el = frac_pT_el.GetXaxis().GetNbins()
nbinsy_el = frac_pT_el.GetYaxis().GetNbins()
wsums_el = []

for binx in range(nbinsx_el):
    wsum=0
    for biny in range(nbinsy_el):
        wsum +=frac_pT_el.GetBinContent(binx,biny)
    wsums_el.append(wsum)
x_pT_el = []
y_frac_el = []
error_el = []
for binx in range(nbinsx_el):
    mean=0
    mean_squared=0
    if wsums_el[binx]==0.0:
        continue
    for biny in range(nbinsy_el):
        bincenter_y = frac_pT_el.GetYaxis().GetBinCenter(biny)
        weight = frac_pT_el.GetBinContent(binx,biny)/wsums_el[binx]
        mean += weight*bincenter_y
        mean_squared +=weight*bincenter_y**2
    y_frac_el.append(mean)
    x_pT_el.append(frac_pT_el.GetXaxis().GetBinCenter(binx))
    error_el.append(math.sqrt(mean_squared-mean**2)/math.sqrt(wsums_el[binx])) #standard error                                                                                        

el_graph = ROOT.TGraphErrors(len(x_pT_el))
el_graph.SetName("el_graph")

assert len(x_pT_el) == len(y_frac_el)
for p in range(len(x_pT_el)):
    el_graph.SetPoint(p, x_pT_el[p], y_frac_el[p])
    el_graph.SetPointError(p, 0, error_el[p])

#------------------------------------------------#
# SAVE HISTOGRAMS
#------------------------------------------------#
output_file = ROOT.TFile("output.root","RECREATE")

h_mu_eta.Write()
h_mu_phi.Write()

cone_rad_mu.SetXTitle("#Delta R")
cone_rad_mu.Write()
    
cone_rad_el.SetXTitle("#Delta R")
cone_rad_el.Write()

pixel_mu.SetXTitle("#eta")
pixel_mu.SetYTitle("#phi")
pixel_mu.Write()

pixel_el.SetXTitle("#eta")
pixel_el.SetYTitle("#phi")
pixel_el.Write()

frac_pT_mu.Write()
frac_pT_el.Write()

mu_graph.GetXaxis().SetTitle("p_{T}")
mu_graph.GetYaxis().SetTitle("Fractional p_{T}")
mu_graph.SetTitle("Fractional p_{T} resolution for muons")
output_file.WriteTObject(mu_graph)

el_graph.GetXaxis().SetTitle("p_{T}")
el_graph.GetYaxis().SetTitle("Fractional p_{T}")
el_graph.SetTitle("Fractional p_{T} resolution for electrons")
el_graph.Write()

output_file.Close()
##################################################
