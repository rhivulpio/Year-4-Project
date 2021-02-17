#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <TCanvas.h>
#include "TF1.h"
#include "TEfficiency.h"
#include <TString.h>

void signal_bkgd(){
    TString signal_file = "signal.root";
    TString bkgd_file = "bkgd.root";

    TFile * signal = TFile::Open(signal_file);
    TFile * bkgd = TFile::Open(bkgd_file);
    TFile * OutputFile = new TFile("signal_bkgd.root", "recreate");
    OutputFile -> cd();
    TCanvas * c1 = new TCanvas("c1"," Efficiency ", 50, 50, 1000, 750);
    gStyle -> SetOptStat(0); //what does this do

    TH1D * h_signal;
    TH1D * h_bkgd;
    TLegend * legend;
    Double_t histmax; //what is this

    
    std::vector<TString> signal_histograms;

    signal_histograms.push_back("Mass Reconstruction/h_Higgs_reco");
    signal_histograms.push_back("Mass Reconstruction/h_Z_reco");
    signal_histograms.push_back("Mass Reconstruction/h_Zstar_reco");

    for(Int_t i = 0; i < signal_histograms.size() ; i++){ //int_t????
        h_signal = (TH1D*) signal -> Get(signal_histograms[i]);
        h_bkgd = (TH1D*) bkgd -> Get(signal_histograms[i]);    
        
        histmax = h_signal -> GetMaximum();
        if(histmax < h_bkgd -> GetMaximum()) histmax = h_bkgd -> GetMaximum();

        h_signal -> SetAxisRange(0, histmax*1.1, "Y");
        h_signal -> SetLineColor(2);

        h_signal -> Draw("hist E2");
        h_bkgd -> Draw("same hist E2");

        legend = new TLegend(0.1, 0.8, 0.25, 0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h_signal, "Signal");
        legend->AddEntry(h_bkgd, "Background");
        legend-> Draw("same"); 
        c1 -> Write(signal_histograms[i]);
    }

    OutputFile -> Close();

}
