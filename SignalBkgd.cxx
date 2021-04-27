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
#include "Utilities.cxx"
#include "GlobalConstants.cxx"

TFile * signal_file;
TFile * bkgd_file;
TFile * bkgd_file_two;
TFile * OutputFile;

void SignalBkgd(){
    TString signal_file_name = "signal.root";
    TString bkgd_file_name = "bkgd.root";
    TString bkgd_file_two_name = "bkgd_two.root";

    signal_file = TFile::Open(signal_file_name);
    bkgd_file = TFile::Open(bkgd_file_name);
    bkgd_file_two = TFile::Open(bkgd_file_two_name);
    OutputFile = new TFile("SignalBkgd.root", "recreate");
    OutputFile -> cd();
    TCanvas * c2 = new TCanvas("c2"," Efficiency ", 50, 50, 1000, 750);
    gStyle -> SetOptStat(0); //what does this do

    TH1D * h_signal;
    TH1D * h_bkgd;
    TH1D * h_bkgd_two;
    TLegend * legend;
    Double_t histmax;
    
    std::vector<TString> histograms;

    histograms.push_back("Mass Reconstruction/h_Higgs_reco");
    histograms.push_back("Mass Reconstruction/h_Z_reco");
    histograms.push_back("Mass Reconstruction/h_Zstar_reco");

    std::vector<TString> smeared_histograms;
    smeared_histograms.push_back("Smearing/h_Higgs_reco_smeared");
    smeared_histograms.push_back("Smearing/h_Z_reco_smeared");
    smeared_histograms.push_back("Smearing/h_Zstar_reco_smeared");

    TH1D * h_signal_data;
    TH1D * h_bkgd_data;

    h_signal_data = (TH1D*) signal_file -> Get(histograms[1]);
    h_bkgd_data = (TH1D*) bkgd_file -> Get(histograms[1]);

    int min = h_signal_data -> FindBin(20);
    int max = h_signal_data -> FindBin(120);
    double Z_signal = h_signal_data -> GetEntries();
    double Z_bkgd = h_bkgd_data -> GetEntries();

    std::cout << "Z signal data = " << Z_signal << " Z background data = " << Z_bkgd << std::endl;

    for(Int_t i = 0; i < histograms.size() ; ++i){
        h_signal = (TH1D*) signal_file -> Get(histograms[i]);
        h_bkgd = (TH1D*) bkgd_file -> Get(histograms[i]);    
        
        histmax = h_signal -> GetMaximum();
        if(histmax < h_bkgd -> GetMaximum()){
            histmax = h_bkgd -> GetMaximum();
        }

        h_signal -> SetAxisRange(0, histmax*1.1, "Y");
        h_signal -> SetLineColor(2);

        h_signal -> Draw("hist E2");
        h_bkgd -> Draw("same hist E2");

        legend = new TLegend(0.1, 0.8, 0.25, 0.9);
        legend->SetHeader("Histogram Markers","C"); // option "C" allows to center the header
        legend->AddEntry(h_signal, "Signal");
        legend->AddEntry(h_bkgd, "Background");
        legend-> Draw("same"); 
        c2 -> Write(histograms[i]);
    }

    //----------------------------------------------------
    //Plotting Significance and Signal to Background Ratio
    //----------------------------------------------------
    Significance_Plots("pT", "GeV", n_pT_cuts, min_pT_cut, step_pT_cut);
    Significance_Plots("Zstar", "GeV", n_Zstar_cuts, min_Zstar_cut, step_Zstar_cut);
    Significance_Plots("Z", "GeV", n_Z_cuts, min_Z_cut, step_Z_cut);

    //hstack(histograms[0]);
    hstack(smeared_histograms[0]);

    OutputFile -> Close();
}

void Significance_Plots(TString property, TString units, int n_cuts, double min_cut, double step_cut){
    std::vector<TString> h_varycuts = Generate_Histogram_List(n_cuts, property);

    TH1D * h_signif = new TH1D("h_signif_"+property, "; " + property + " Cut [" + units + "]; Significance;", n_cuts, min_cut, min_cut + n_cuts*step_cut);
    TH1D * h_ratio = new TH1D("h_ratio_"+property, "; " + property + " Cut [" + units + "]; Signal to Background Ratio;", n_cuts, min_cut, min_cut + n_cuts*step_cut);

    Int_t higgs_min;
    Int_t higgs_max;
    double signal_data;
    double bkgd_data;
    double significance;
    double ratio;
    int bin_pt;
    TH1D * h_signal_cuts;
    TH1D * h_bkgd_cuts;  

    std::cout << h_varycuts[0] << std::endl;

    for(int n = 0; n < n_cuts; ++n){
        h_signal_cuts = (TH1D*) signal_file -> Get(h_varycuts[n]);
        h_bkgd_cuts = (TH1D*) bkgd_file -> Get(h_varycuts[n]);

        higgs_min = h_signal_cuts -> FindBin(120);
        higgs_max = h_signal_cuts -> FindBin(130);
        signal_data = h_signal_cuts -> Integral(higgs_min, higgs_max);
        bkgd_data = h_bkgd_cuts -> Integral(higgs_min, higgs_max);
        significance = signal_data/TMath::Sqrt(bkgd_data);
        ratio = signal_data/bkgd_data;

        // std::cout << "significance = " << significance << std::endl;
        // std::cout << "signal to background ratio = " << ratio << std::endl;

        std::vector<double> cut_values = Define_Cut_Values(n_cuts, min_cut, step_cut);

        bin_pt = h_signif -> FindBin(cut_values.at(n));

        bool infinity_sig = isinf(significance);
        bool nan_sig = isnan(significance);
        bool infinity_ratio = isinf(ratio);
        bool nan_ratio = isnan(ratio);

        if(!infinity_sig && !nan_sig && !infinity_ratio && !nan_ratio){
            h_signif -> SetBinContent(bin_pt, significance);
            h_ratio -> SetBinContent(bin_pt, ratio);
        }
    }

    OutputFile->cd();
    h_signif->Write();
    h_ratio->Write();

}  

TCanvas * hstack(TString histogram) {
    THStack * hs = new THStack("hs", "Stacked Mass Reconstruction");
    
    TH1D * h2 = (TH1D*) bkgd_file -> Get(histogram); 
    h2->SetFillColor(29);
    h2->SetLineColor(29);
    //h2->SetFillStyle(3001);
    hs->Add(h2);

    TH1D * h3 = (TH1D*) bkgd_file_two -> Get(histogram); 
    h3->SetFillColor(46);
    h3->SetLineColor(46);
    //h3->SetFillStyle(3001);
    hs->Add(h3);
    
    TH1D * h1 = (TH1D*) signal_file -> Get(histogram);
    h1->SetFillColor(38);
    h1->SetLineColor(38);
    //h1->SetFillStyle(3001);
    // h1->GetXaxis()->SetRangeUser(80, 170);
    hs->Add(h1);

    

    TCanvas * cst = new TCanvas("cst", "Stacked Mass Reconstruction", 50, 50, 1000, 750);
    // h1 -> Draw();
    // h2 -> Draw();
    hs->Draw("HIST");

    hs->GetXaxis()->SetTitle("m_{4l} (GeV)");
    hs->GetYaxis()->SetTitle("Number of Events");
    //hs->GetXaxis()->SetRangeUser(0, 200);

    TLegend * legend2;
    legend2 = new TLegend(0.1, 0.7, 0.3, 0.9);
    // legend->SetHeader("Particle", "C");
    legend2->AddEntry(h1, "Signal");
    legend2->AddEntry(h2, "ZZ* -> 4l Background");
    legend2->AddEntry(h3, "Z -> 4l Background");
    legend2->Draw("same");

    hs->Write();

    return cst;
}