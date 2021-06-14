
import sys,os,re,time,math
import array
from time import sleep
import threading
import multiprocessing
#import random
import math
from math import sqrt
import yaml

import ROOT
from ROOT import TGraph,TGraphErrors,TH1F,TCanvas,TFile,TLegend,TLatex

PP_MT_VAL = [1.10765, 1.16831, 1.22838, 1.31558, 1.46282, 1.68716, 2.21156]
PL_MT_VAL = [1.2124, 1.2896, 1.376, 1.5407, 1.756, 2.2594]

def DrawErrorBand(graph):
    isErrorBand = graph.GetErrorYhigh(0) != -1 and graph.GetErrorYlow(0) != -1
    npoints     = graph.GetN()

    if not isErrorBand:
        graph.Draw("l same")
        return

    # Declare individual TGraph objects used in drawing error band
    central, min, max = ROOT.TGraph(), ROOT.TGraph(), ROOT.TGraph()
    shapes = []
    for i in range((npoints-1)*4):
        shapes.append(ROOT.TGraph())

    # Set ownership of TGraph objects
    ROOT.SetOwnership(central, False)
    ROOT.SetOwnership(    min, False)
    ROOT.SetOwnership(    max, False)
    for shape in shapes:
        ROOT.SetOwnership(shape, False)

    # Get data points from TGraphAsymmErrors
    x, y, ymin, ymax = [], [], [], []
    for i in range(npoints):
        tmpX = array.array('d', [0])
        tmpY = array.array('d', [0])
        graph.GetPoint(i, tmpX, tmpY)
        x.append(tmpX[0])
        y.append(tmpY[0])
        ymin.append(tmpY[0] - graph.GetErrorYlow(i))
        ymax.append(tmpY[0] + graph.GetErrorYhigh(i))

    # Fill central, min and max graphs
    for i in range(npoints):
        central.SetPoint(i, x[i], y[i])
        min.SetPoint(i, x[i], ymin[i])
        max.SetPoint(i, x[i], ymax[i])

    # Fill shapes which will be shaded to create the error band
    for i in range(npoints-1):
        for version in range(4):
            shapes[i+(npoints-1)*version].SetPoint((version+0)%4, x[i],   ymax[i])
            shapes[i+(npoints-1)*version].SetPoint((version+1)%4, x[i+1], ymax[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+2)%4, x[i+1], ymin[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+3)%4, x[i],   ymin[i])

    # Set attributes to those of input graph
    central.SetLineColor(graph.GetLineColor())
    central.SetLineStyle(graph.GetLineStyle())
    central.SetLineWidth(graph.GetLineWidth())
    min.SetLineColor(graph.GetLineColor())
    min.SetLineStyle(graph.GetLineStyle())
    max.SetLineColor(graph.GetLineColor())
    max.SetLineStyle(graph.GetLineStyle())
    for shape in shapes:
        shape.SetFillColor(graph.GetFillColor())
        shape.SetFillStyle(graph.GetFillStyle())

    # Draw
    for shape in shapes:
        shape.Draw("f same")
    min.Draw("l same")
    max.Draw("l same")
    central.Draw("l same")
    ROOT.gPad.RedrawAxis()
    sleep(10)

def ErrorBandInMt(InputFilePath,mT,DATA,SourceVar,SmearStrategy,system):
    OutputPDF = InputFilePath+'cFinal_'+SourceVar+'_SS'+str(SmearStrategy)+'_mT'+str(mT)+'.pdf'
    OutputROOT = InputFilePath+'QA_'+SourceVar+'_SS'+str(SmearStrategy)+'_mT'+str(mT)+'.root'
    mean_rad=0
    stdv_rad=0
    stat_rad=0
    n_rad=0;

    kstar_data_low = []
    kstar_data_up = []
    kstar_data = []
    mean_data = []
    stdv_data = []
    stat_data = []#evaluated from the mean error of the bin
    n_data = []
    kstar_fit = []
    mean_fit = []
    stdv_fit = []
    n_fit = []

    hQA = TH1F("hQA_mT"+str(mT),"hQA_mT"+str(mT),512,0.5,2.0)

    for data in DATA:
        InputFileName = InputFilePath+'fOut_'+SourceVar+'_SS'+str(SmearStrategy)+'_D'+str(data)+'_mT'+str(mT)+'.root'
        infile = ROOT.TFile.Open(InputFileName)
        TreeName = 'ppTree'
        if "pLambda" in system:
            TreeName='pLTree'
        ppTree = infile.Get(TreeName)

        for event in ppTree :
            #gData = ROOT.TGraph()
            hData = event.hData
            gFit = event.gFit
            #hData.Draw()
            #sleep(11)
            #IF YOU HAVE SOME CONDITIONS TO REJECT, THIS WOULD BE THE PLACE
            if event.BaselineVar>1:
                continue
            if system=="pLambda_LO" and (event.pL_pot_var!=-11600):
                continue
            if system=="pLambda_NLO" and (event.pL_pot_var!=11600):
                continue
            #if event.FemtoRegion!=376:
            #    continue
            #----------------------------------------------------
            mean_rad += event.rad
            hQA.Fill(event.rad)
            stdv_rad += event.rad**2
            stat_rad += event.raderr
            n_rad += 1
            for point in range(0,gFit.GetN()):
                kstar = array.array('d', [0])
                valY = array.array('d', [0])
                gFit.GetPoint(point, kstar, valY)
                if len(mean_fit)<=point:
                    mean_fit.append(valY[0])
                    stdv_fit.append(valY[0]*valY[0])
                    kstar_fit.append(kstar[0])
                    n_fit.append(1)
                else:
                    mean_fit[point] = float(n_fit[point])/float(n_fit[point]+1)*(mean_fit[point]+valY[0]/n_fit[point])
                    stdv_fit[point] = float(n_fit[point])/float(n_fit[point]+1)*(stdv_fit[point]+valY[0]*valY[0]/n_fit[point])
                    kstar_fit[point] = kstar[0]
                    n_fit[point] += 1
            for bin in range(0,hData.GetNbinsX()):
                kstar = hData.GetBinCenter(bin+1)
                valY = hData.GetBinContent(bin+1)
                errY = hData.GetBinError(bin+1)
                if len(mean_data)<=bin:
                    mean_data.append(valY)
                    stdv_data.append(valY*valY)
                    kstar_data.append(kstar)
                    kstar_data_low.append(hData.GetXaxis().GetBinLowEdge(bin+1))
                    kstar_data_up.append(hData.GetXaxis().GetBinUpEdge(bin+1))
                    stat_data.append(errY)
                    n_data.append(1)
                else:
                    mean_data[bin] = float(n_data[bin])/float(n_data[bin]+1)*(mean_data[bin]+valY/n_data[bin])
                    stdv_data[bin] = float(n_data[bin])/float(n_data[bin]+1)*(stdv_data[bin]+valY*valY/n_data[bin])
                    kstar_data[bin] = kstar
                    #kstar_data_low,kstar_data_up
                    stat_data[bin] = float(n_data[bin])/float(n_data[bin]+1)*(stat_data[bin]+errY/n_data[bin])
                    n_data[bin] += 1

    mean_rad /= float(n_rad)
    stdv_rad /= float(n_rad)
    stdv_rad = sqrt(stdv_rad-mean_rad**2)
    stat_rad /= float(n_rad)
    #print "rad = {:.4f} (stat) {:.4f} (syst) {:.4f} fm".format(mean_rad,stat_rad,stdv_rad)

    hFinalDataStat = TH1F("hFinalDataStat","hFinalDataStat",len(n_data),kstar_data_low[0],kstar_data_up[-1])
    hFinalDataStat.SetStats(0);
    hFinalDataStat.SetTitle("");
    hFinalDataStat.GetXaxis().SetTitle("k* (MeV)");
    hFinalDataStat.GetXaxis().SetTitleSize(0.06);
    hFinalDataStat.GetXaxis().SetLabelSize(0.06);
    hFinalDataStat.GetXaxis().CenterTitle();
    hFinalDataStat.GetXaxis().SetTitleOffset(1.3);
    hFinalDataStat.GetXaxis().SetLabelOffset(0.02);
    hFinalDataStat.GetYaxis().SetTitle("C(k*)");
    hFinalDataStat.GetYaxis().SetTitleSize(0.06);
    hFinalDataStat.GetYaxis().SetLabelSize(0.06);
    hFinalDataStat.GetYaxis().CenterTitle();
    hFinalDataStat.GetYaxis().SetTitleOffset(1.10);
    if "pLambda" in system:
        hFinalDataStat.GetYaxis().SetRangeUser(0.9, 2.1);
    else:
        hFinalDataStat.GetYaxis().SetRangeUser(0.9, 2.9);
    hFinalDataStat.GetXaxis().SetRangeUser(0, 200);
    hFinalDataStat.SetLineColor(ROOT.kBlack)
    hFinalDataStat.SetLineWidth(2)
    for bin in range(0,len(n_data)):
        if stdv_data[bin]<=mean_data[bin]**2:
            stdv_data[bin]=0
        else:
            stdv_data[bin] = sqrt(stdv_data[bin]-mean_data[bin]**2)
        #print (stdv_data[bin])
        hFinalDataStat.SetBinContent(bin+1,mean_data[bin])
        hFinalDataStat.SetBinError(bin+1,stat_data[bin])

    gFinalFit = TGraphErrors()
    gFinalFit.SetName("gFinalFit");
    gFinalFit.Set(len(n_fit));
    gFinalFit.SetFillColorAlpha(ROOT.kBlue+1,0.5);
    gFinalFit.SetLineColor(ROOT.kBlue+1);
    gFinalFit.SetLineWidth(1);
    for point in range(0,len(n_fit)):
        stdv_fit[point] = sqrt(stdv_fit[point]-mean_fit[point]**2)
        gFinalFit.SetPoint(point,kstar_fit[point],mean_fit[point])
        gFinalFit.SetPointError(point,2,stdv_fit[point])

    fQA = TFile (OutputROOT,"recreate")
    hQA.Write()
    #print (str(stdv_rad)+' vs '+str(hQA.GetStdDev()))
    #cFinal = TCanvas("cFinal", "cFinal", 1);
    #cFinal.cd(0)
    #cFinal.SetCanvasSize(1280, 720)
    #cFinal.SetMargin(0.15,0.05,0.2,0.05);
    #hFinalDataStat.Draw()
    #gFinalFit.Draw("same 3l")
    #cFinal.SaveAs(OutputPDF)
    return mean_rad,stat_rad,stdv_rad
    #input("Press Enter to continue...")
    #for point in range(0,len(kstar_fit)):
    #    stdv_fit[point] = sqrt(stdv_fit[point]-mean_fit[point]**2)
    #    print 'Point at '+"{:.0f}".format(kstar_fit[point])+' is '+"{:.3f}".format(mean_fit[point])+" +/- {:.3f}".format(stdv_fit[point])


def mT_dependance(InputFilePath,src,MT_BINS,MT_VAL,DATA,system):
    gMtStat = TGraphErrors()
    gMtStat.SetName("gMtStat")

    gMtSyst = TGraphErrors()
    gMtSyst.SetName("gMtSyst")

    for mT in MT_BINS:
        rad,raderr,radsys = ErrorBandInMt(InputFilePath,mT,DATA,src,1,system)
        gMtStat.SetPoint(mT,MT_VAL[mT],rad)
        gMtStat.SetPointError(mT,0.0,raderr)
        gMtSyst.SetPoint(mT,MT_VAL[mT],rad)
        gMtSyst.SetPointError(mT,0.0,radsys)
    return gMtStat, gMtSyst

#print 'Hello'
#def ErrorBandInMt(InputFilePath,mT,DATA,SourceVar,SmearStrategy):
#ErrorBandInMt('/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/Systematics/',0,[0],'Gauss',1)
#ErrorBandInMt('/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/Systematics/',0,[0,1,2,3,4],'McGauss_ResoTM',1)

def Make_pp_plots():
    PP_PATH = '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/Systematics/'
    OutFileName = PP_PATH+"gMt_pp.root"
    outputFile = TFile(OutFileName,"recreate")
    gMtStatGauss_pp, gMtSystGauss_pp = mT_dependance(PP_PATH,"Gauss", range(0,7), PP_MT_VAL, range(0,27), "pp")
    gMtStatCore_pp, gMtSystCore_pp = mT_dependance(PP_PATH,"McGauss_ResoTM", range(0,7), PP_MT_VAL, range(0,27), "pp")
    outputFile.cd()
    gMtStatGauss_pp.SetName("gMtStatGauss_pp")
    gMtStatGauss_pp.Write()
    gMtSystGauss_pp.SetName("gMtSystGauss_pp")
    gMtSystGauss_pp.Write()
    gMtStatCore_pp.SetName("gMtStatCore_pp")
    gMtStatCore_pp.Write()
    gMtSystCore_pp.SetName("gMtSystCore_pp")
    gMtSystCore_pp.Write()

    return OutFileName

def Make_pL_plots():
    PL_PATH = '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pL/Systematics/'
    OutFileName = PL_PATH+"gMt_pL.root"
    outputFile = TFile(OutFileName,"recreate")
    gMtStatGauss_pL_LO, gMtSystGauss_pL_LO = mT_dependance(PL_PATH,"Gauss", range(0,6), PL_MT_VAL, range(0,43), "pLambda_LO")
    gMtStatCore_pL_LO, gMtSystCore_pL_LO = mT_dependance(PL_PATH,"McGauss_ResoTM", range(0,6), PL_MT_VAL, range(0,43), "pLambda_LO")
    outputFile.cd()
    gMtStatGauss_pL_LO.SetName("gMtStatGauss_pL_LO")
    gMtStatGauss_pL_LO.Write()
    gMtSystGauss_pL_LO.SetName("gMtSystGauss_pL_LO")
    gMtSystGauss_pL_LO.Write()
    gMtStatCore_pL_LO.SetName("gMtStatCore_pL_LO")
    gMtStatCore_pL_LO.Write()
    gMtSystCore_pL_LO.SetName("gMtSystCore_pL_LO")
    gMtSystCore_pL_LO.Write()

    gMtStatGauss_pL_NLO, gMtSystGauss_pL_NLO = mT_dependance(PL_PATH,"Gauss", range(0,6), PL_MT_VAL, range(0,43), "pLambda_NLO")
    gMtStatCore_pL_NLO, gMtSystCore_pL_NLO = mT_dependance(PL_PATH,"McGauss_ResoTM", range(0,6), PL_MT_VAL, range(0,43), "pLambda_NLO")
    outputFile.cd()
    gMtStatGauss_pL_NLO.SetName("gMtStatGauss_pL_NLO")
    gMtStatGauss_pL_NLO.Write()
    gMtSystGauss_pL_NLO.SetName("gMtSystGauss_pL_NLO")
    gMtSystGauss_pL_NLO.Write()
    gMtStatCore_pL_NLO.SetName("gMtStatCore_pL_NLO")
    gMtStatCore_pL_NLO.Write()
    gMtSystCore_pL_NLO.SetName("gMtSystCore_pL_NLO")
    gMtSystCore_pL_NLO.Write()

    return OutFileName

#3 plots, one for Gaus and one for Core
#one with both results combined together (also adding NLO and LO), where we plot as a gray
#background the published results for comparison
def Make_mT_scaling_plots(ppFileName,pLFileName,HepDataFile,OutputFolder):
    dAxis = [
            {"histo":TH1F("hAxis_core","hAxis_core",128,0.95,2.45),"rMin":0.67,"rMax":1.43,"yaxis":"#it{r}_{core} (fm)"},
            {"histo":TH1F("hAxis_gaus","hAxis_gaus",128,0.95,2.45),"rMin":0.87,"rMax":1.65,"yaxis":"#it{r}_{0} (fm)"},
            {"histo":TH1F("hAxis_both","hAxis_both",128,0.95,2.45),"rMin":0.67,"rMax":1.65,"yaxis":"#it{r}_{core},#it{r}_{0} (fm)"}
            ]
    #dAxis[0]["histo"].SetTitle("; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); #it{r}_{core} (fm)")
    for axis in dAxis:
        axis["histo"].SetTitle("; #LT#it{m}_{T}#GT  (GeV/#it{c}^{2}); "+axis["yaxis"])
        axis["histo"].SetStats(False)
        axis["histo"].GetXaxis().SetTitleSize(0.06);
        axis["histo"].GetYaxis().SetTitleSize(0.06);
        axis["histo"].GetXaxis().SetTitleOffset(1.15);
        axis["histo"].GetYaxis().SetTitleOffset(1.1);
        axis["histo"].GetXaxis().SetLabelSize(0.06);
        axis["histo"].GetYaxis().SetLabelSize(0.06);
        axis["histo"].GetXaxis().SetLabelOffset(0.02);
        axis["histo"].GetYaxis().SetLabelOffset(0.02);
        axis["histo"].GetXaxis().SetRangeUser(0.95, 2.7);
        axis["histo"].GetYaxis().SetRangeUser(axis["rMin"], axis["rMax"]);
        #for bin in range(0,axis["histo"].GetNbinsX()):
        #    axis["histo"].SetBinContent(bin+1,1)

    tFile = []
    tFile.append(TFile(ppFileName,"read"))#0
    tFile.append(TFile(pLFileName,"read"))#1
    tFile.append(TFile(HepDataFile,"read"))#2

    #3225
    dResult = [
    {"system":"pp","description":"p#minusp (AV18)","version":"PLB","source":"Gauss","histoStat":tFile[2].geHEP_pp_Gaus_Stat,"histoSyst":tFile[2].geHEP_pp_Gaus_Syst,
    "color_main":ROOT.kBlue+2,"color_combi":ROOT.kBlue+2,"color_shady":ROOT.kGray+2,"marker_style":34,"marker_size":2.0,"fill_style":3244},
    {"system":"pp","description":"p#minusp (AV18)","version":"PLB","source":"McGauss_ResoTM","histoStat":tFile[2].geHEP_pp_Core_Stat,"histoSyst":tFile[2].geHEP_pp_Core_Syst,
    "color_main":ROOT.kBlue+2,"color_combi":ROOT.kBlue+2,"color_shady":ROOT.kGray+2,"marker_style":34,"marker_size":2.0,"fill_style":3244},
    {"system":"pp","description":"p#minusp (AV18)","version":"Ana21v1","source":"Gauss","histoStat":tFile[0].gMtStatGauss_pp,"histoSyst":tFile[0].gMtSystGauss_pp,
    "color_main":ROOT.kBlue+2,"color_combi":ROOT.kBlue+2,"color_shady":ROOT.kGray+2,"marker_style":34,"marker_size":2.0,"fill_style":3244},
    {"system":"pp","description":"p#minusp (AV18)","version":"Ana21v1","source":"McGauss_ResoTM","histoStat":tFile[0].gMtStatCore_pp,"histoSyst":tFile[0].gMtSystCore_pp,
    "color_main":ROOT.kBlue+2,"color_combi":ROOT.kBlue+2,"color_shady":ROOT.kGray+2,"marker_style":34,"marker_size":2.0,"fill_style":3244},

    {"system":"pLambdaLO","description":"p#minus#Lambda (LO)","version":"PLB","source":"Gauss","histoStat":tFile[2].geHEP_pL_LO_Gaus_Stat,"histoSyst":tFile[2].geHEP_pL_LO_Gaus_Syst,
    "color_main":ROOT.kGreen+3,"color_combi":ROOT.kGreen+3,"color_shady":ROOT.kGray+2,"marker_style":49,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaLO","description":"p#minus#Lambda (LO)","version":"PLB","source":"McGauss_ResoTM","histoStat":tFile[2].geHEP_pL_LO_Core_Stat,"histoSyst":tFile[2].geHEP_pL_LO_Core_Syst,
    "color_main":ROOT.kGreen+3,"color_combi":ROOT.kGreen+3,"color_shady":ROOT.kGray+2,"marker_style":49,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaLO","description":"p#minus#Lambda (LO)","version":"Ana21v1","source":"Gauss","histoStat":tFile[1].gMtStatGauss_pL_LO,"histoSyst":tFile[1].gMtSystGauss_pL_LO,
    "color_main":ROOT.kGreen+3,"color_combi":ROOT.kGreen+3,"color_shady":ROOT.kGray+2,"marker_style":49,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaLO","description":"p#minus#Lambda (LO)","version":"Ana21v1","source":"McGauss_ResoTM","histoStat":tFile[1].gMtStatCore_pL_LO,"histoSyst":tFile[1].gMtSystCore_pL_LO,
    "color_main":ROOT.kGreen+3,"color_combi":ROOT.kGreen+3,"color_shady":ROOT.kGray+2,"marker_style":49,"marker_size":2.0,"fill_style":3244},

    {"system":"pLambdaNLO","description":"p#minus#Lambda (NLO)","version":"PLB","source":"Gauss","histoStat":tFile[2].geHEP_pL_NLO_Gaus_Stat,"histoSyst":tFile[2].geHEP_pL_NLO_Gaus_Syst,
    "color_main":ROOT.kRed+1,"color_combi":ROOT.kRed+1,"color_shady":ROOT.kGray+2,"marker_style":41,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaNLO","description":"p#minus#Lambda (NLO)","version":"PLB","source":"McGauss_ResoTM","histoStat":tFile[2].geHEP_pL_NLO_Core_Stat,"histoSyst":tFile[2].geHEP_pL_NLO_Core_Syst,
    "color_main":ROOT.kRed+1,"color_combi":ROOT.kRed+1,"color_shady":ROOT.kGray+2,"marker_style":41,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaNLO","description":"p#minus#Lambda (NLO)","version":"Ana21v1","source":"Gauss","histoStat":tFile[1].gMtStatGauss_pL_NLO,"histoSyst":tFile[1].gMtSystGauss_pL_NLO,
    "color_main":ROOT.kRed+1,"color_combi":ROOT.kRed+1,"color_shady":ROOT.kGray+2,"marker_style":41,"marker_size":2.0,"fill_style":3244},
    {"system":"pLambdaNLO","description":"p#minus#Lambda (NLO)","version":"Ana21v1","source":"McGauss_ResoTM","histoStat":tFile[1].gMtStatCore_pL_NLO,"histoSyst":tFile[1].gMtSystCore_pL_NLO,
    "color_main":ROOT.kRed+1,"color_combi":ROOT.kRed+1,"color_shady":ROOT.kGray+2,"marker_style":41,"marker_size":2.0,"fill_style":3244}
    ]
#lbrt
    Canvas = [
    {"caption":"Gaussian Source (PLB)","canvas":TCanvas("c_PLB_Gauss","c_PLB_Gauss",1200,800),"name":"c_PLB_Gauss","legend":TLegend(0.18, 0.23, 0.46, 0.49),
    "MainHisto":
    [{"system":"pp","version":"PLB","source":"Gauss"}, {"system":"pLambdaLO","version":"PLB","source":"Gauss"}, {"system":"pLambdaNLO","version":"PLB","source":"Gauss"}],
    "ShadyHisto":
    []},
    {"caption":"Gaussian Source (NEW)","canvas":TCanvas("c_Ana21v1_Gauss","c_Ana21v1_Gauss",1200,800),"name":"c_Ana21v1_Gauss","legend":TLegend(0.18, 0.23, 0.46, 0.49),
    "MainHisto":
    [{"system":"pp","version":"Ana21v1","source":"Gauss"}, {"system":"pLambdaLO","version":"Ana21v1","source":"Gauss"}, {"system":"pLambdaNLO","version":"Ana21v1","source":"Gauss"}],
    "ShadyHisto":
    []},

    {"caption":"Gaussian + Resonances (PLB)","canvas":TCanvas("c_PLB_Core","c_PLB_Core",1200,800),"name":"c_PLB_Core","legend":TLegend(0.18, 0.23, 0.46, 0.49),
    "MainHisto":
    [{"system":"pp","version":"PLB","source":"McGauss_ResoTM"}, {"system":"pLambdaLO","version":"PLB","source":"McGauss_ResoTM"}, {"system":"pLambdaNLO","version":"PLB","source":"McGauss_ResoTM"}],
    "ShadyHisto":
    []},
    {"caption":"Gaussian + Resonances (NEW)","canvas":TCanvas("c_Ana21v1_Core","c_Ana21v1_Core",1200,800),"name":"c_Ana21v1_Core","legend":TLegend(0.18, 0.23, 0.46, 0.49),
    "MainHisto":
    [{"system":"pp","version":"Ana21v1","source":"McGauss_ResoTM"}, {"system":"pLambdaLO","version":"Ana21v1","source":"McGauss_ResoTM"}, {"system":"pLambdaNLO","version":"Ana21v1","source":"McGauss_ResoTM"}],
    "ShadyHisto":
    []}
    ]

    final_pdf = []

    for can in Canvas:
        can["canvas"].SetMargin(0.15,0.05,0.2,0.05)
        can["legend"].SetFillStyle(0)
        can["legend"].SetTextFont(43)
        can["legend"].SetTextSizePixels(40)
        can["legend"].SetBorderSize(0)
        NumGraphs = len(can["MainHisto"])
        histoStat = []
        histoSyst = []
        system = []
        version = []
        source = []
        SameSourceType = True
        for i in range(0,NumGraphs):
            system.append(can["MainHisto"][i]["system"])
            version.append(can["MainHisto"][i]["version"])
            source.append(can["MainHisto"][i]["source"])
            for result in dResult:
                for hist in ["histoStat","histoSyst"]:
                    result[hist].SetLineColor(result["color_main"])
                    result[hist].SetLineWidth(2)
                    result[hist].SetFillColor(result["color_main"])
                    result[hist].SetMarkerColor(result["color_main"])
                    result[hist].SetFillColorAlpha(result["color_main"], 0.5)
                    #result[hist].SetFillStyle(result["fill_style"])
                    result[hist].SetMarkerStyle(result["marker_style"])
                    result[hist].SetMarkerSize(result["marker_size"])
                if result["system"]==system[i] and result["version"]==version[i] and result["source"]==source[i]:
                    histoStat.append(result["histoStat"])
                    histoSyst.append(result["histoSyst"])
                    can["legend"].AddEntry(result["histoStat"],result["description"])
            if i>0 and source[i-1]!=source[i]:
                SameSourceType = False
        #by now we have found all histograms we want to plot
        if SameSourceType==False:
            axisID = 2
        elif source[0]=="Gauss":
            axisID = 1
        else:
            axisID = 0
        can["canvas"].cd()
        dAxis[axisID]["histo"].Draw("axis")
        for hst in histoSyst:
            hst.Draw("e3 same")
        for hst in histoStat:
            hst.Draw("pez same")
        can["legend"].Draw("same")

        BeamText = TLatex();
        BeamText.SetTextFont(43)
        BeamText.SetTextSize(36)
        BeamText.SetNDC(True)
        BeamText.DrawLatex(0.48, 0.86,"ALICE pp #sqrt{#it{s}} = 13 TeV")
        BeamText.DrawLatex(0.48,0.79,
          "High-mult. (0#kern[-0.65]{ }#minus#kern[-0.65]{ }0.17#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)")
        BeamText.DrawLatex(0.48, 0.72, can["caption"])

        CanvasBaseName = OutputFolder+can["name"]
        can["canvas"].SaveAs(CanvasBaseName+".pdf")
        os.system('convert -density 300 -trim '+CanvasBaseName+'.pdf -quality 100 '+CanvasBaseName+'.png')
        final_pdf.append(CanvasBaseName+".pdf")

    OSCOMMAND = 'pdfunite'
    for str in final_pdf:
        OSCOMMAND += ' '
        OSCOMMAND += str
    OSCOMMAND += ' '
    OSCOMMAND += OutputFolder
    OSCOMMAND += 'SourceSize_vs_mT.pdf'
    os.system(OSCOMMAND)
    #outFile = TFile(OutputFolder+"")


    #dResult[0]["histo"].Draw()

    #dResult =   [
    #            {"system":"p#minusp","version":..dimi or bernie,"source":"Gauss","histo":,"color":..,"marker_style":,"marker_size","fill_style":}
    #            ]

def Make_PLB_Resuls():
    FolderPath = "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/PublishedResults/"
    #OutputFileName = FolderPath+"PLBdata.root"
    geHEP = [
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pp_Gaus","yaml":FolderPath+"pp_Gauss.yaml","yamlID":0},
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pp_Core","yaml":FolderPath+"pp_Core.yaml","yamlID":0},
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pL_LO_Gaus","yaml":FolderPath+"pL_Gauss.yaml","yamlID":1},
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pL_LO_Core","yaml":FolderPath+"pL_Core.yaml","yamlID":1},
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pL_NLO_Gaus","yaml":FolderPath+"pL_Gauss.yaml","yamlID":0},
    {"graphStat":TGraphErrors(),"graphSyst":TGraphErrors(),"name":"geHEP_pL_NLO_Core","yaml":FolderPath+"pL_Core.yaml","yamlID":0}
    ]

    #with open('/home/dimihayl/Software/LocalFemto/Output/SourceStudies/PublishedResults/pL_Gauss.yaml') as f:
    #    data = yaml.load(f, Loader=yaml.FullLoader)
    #NumMtBins = len(data["independent_variables"][0]["values"])
    #for bin in range(0,NumMtBins):
    #    mtVal = data["independent_variables"][0]["values"][bin]["value"]
    #    binVal = data["dependent_variables"][1]["values"][bin]["value"]
    #    binSyst = data["dependent_variables"][1]["values"][bin]["errors"][0]["symerror"]
    #    binStat = data["dependent_variables"][1]["values"][bin]["errors"][1]["symerror"]
    #    print (binVal)
        #geHEP_pp_Gaus_Syst.SetPoint(bin,mtVal,binVal)
        #geHEP_pp_Gaus_Syst.SetPointError(bin,0,binSyst)
        #geHEP_pp_Gaus_Stat.SetPoint(bin,mtVal,binVal)
        #geHEP_pp_Gaus_Stat.SetPointError(bin,0,binStat)

    fOutput = TFile(FolderPath+"PLBdata.root","recreate")

    for system in geHEP:
        with open(system["yaml"]) as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        NumMtBins = len(data["independent_variables"][0]["values"])
        for bin in range(0,NumMtBins):
            mtVal = data["independent_variables"][0]["values"][bin]["value"]
            binVal = data["dependent_variables"][system["yamlID"]]["values"][bin]["value"]
            binSyst = data["dependent_variables"][system["yamlID"]]["values"][bin]["errors"][0]["symerror"]
            binStat = data["dependent_variables"][system["yamlID"]]["values"][bin]["errors"][1]["symerror"]
            system["graphStat"].SetName(system["name"]+"_Stat")
            system["graphStat"].SetPoint(bin,mtVal,binVal)
            system["graphStat"].SetPointError(bin,0,binStat)
            system["graphSyst"].SetName(system["name"]+"_Syst")
            system["graphSyst"].SetPoint(bin,mtVal,binVal)
            system["graphSyst"].SetPointError(bin,0,binSyst)
        fOutput.cd()
        system["graphStat"].Write()
        system["graphSyst"].Write()
    #print (len(geHEP))


    #mtval
    #for value in data["independent_variables"][0]["values"]:
    #    print (value["value"])

    #print (data["dependent_variables"][0]["values"][0]["errors"][0]["symerror"])
    #PLB_pp = [
    #{"mT":1.1077,"rad":1.3717,"syst":0.025241,"stat":0.0079898}
    #]
    #geHEP_pp_Gaus_Syst = TGraphErrors()
    #geHEP_pp_Gaus_Syst.SetPoint(0,1.1077,1.3717)
    #geHEP_pp_Gaus_Syst.SetPointError(0,0,0.025241)
    #geHEP_pp_Gaus_Syst.SetPoint(1,1.1683,1.3244)
    #geHEP_pp_Gaus_Syst.SetPointError(1,0,0.022749)
    #geHEP_pp_Gaus_Syst.SetPoint(2,1.1683,1.3244)
    #geHEP_pp_Gaus_Syst.SetPointError(2,0,0.022749)



#ppFileName = Make_pp_plots()
#pLFileName = Make_pL_plots()
#somelist = {'raincoat': 1, 'coins': 5, 'books': 23, 'spectacles': 2, 'chairs': 12, 'pens': 6}
Make_mT_scaling_plots(  '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root',
                        '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root',
                        '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/PublishedResults/PLBdata.root',
                        '/home/dimihayl/Software/LocalFemto/Output/SourceStudies/')
#Make_PLB_Resuls()
