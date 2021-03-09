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
#include "utilities.cxx"

TFile * signal_file;
TFile * bkgd_file;
TFile * OutputFile;

void signal_bkgd(){
    TString signal_file_name = "signal.root";
    TString bkgd_file_name = "bkgd.root";

    signal_file = TFile::Open(signal_file_name);
    bkgd_file = TFile::Open(bkgd_file_name);
    OutputFile = new TFile("signal_bkgd.root", "recreate");
    OutputFile -> cd();
    TCanvas * c1 = new TCanvas("c1"," Efficiency ", 50, 50, 1000, 750);
    gStyle -> SetOptStat(0); //what does this do

    TH1D * h_signal;
    TH1D * h_bkgd;
    TLegend * legend;
    Double_t histmax;

    
    std::vector<TString> histograms;
    //std::vector<TString> h_cuts;

    histograms.push_back("Mass Reconstruction/h_Higgs_reco");
    histograms.push_back("Mass Reconstruction/h_Z_reco");
    histograms.push_back("Mass Reconstruction/h_Zstar_reco");

    //h_cuts.push_back("Mass Reco Cuts/h_Zstar_cuts");

    for(Int_t i = 0; i < histograms.size() ; ++i){
        h_signal = (TH1D*) signal_file -> Get(histograms[i]);
        h_bkgd = (TH1D*) bkgd_file -> Get(histograms[i]);    
        
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
        c1 -> Write(histograms[i]);
    }

    //------------------------
    //calculating significance
    //------------------------
    int n_pT_cuts =40;
    double min_cut = 10.0; //GeV
    double step_cut = 2.0; //GeV
    //std::vector<double> pT_cut_values = Define_Cut_Values(20, 10.0, 2.0);
    //std::vector<TString> h_pT_cuts = Generate_Histogram_List(n_pT_cuts, "pT");
    

    // TH1D * h_signal_cuts;
    // TH1D * h_bkgd_cuts;

    // TH1D * h_signif;
    // TH1D * h_ratio;

    Significance_Plots("pT", "GeV", 40, 10.0, 2.0);
    Significance_Plots("Zstar", "GeV", 40, 10.0, 5.0);

    // Make_Histograms(h_signif_pT, h_ratio_pT, "pT", "GeV", 20, 10.0, 2.0);
    // output = Calculate_Parameters(h_signal_pT, h_bkgd_pT, h_pT_cuts, 20);
    // significance_pT = output[0];
    // ratio_pT = output[1];

    // h_signif_pT = make_histograms_output[0];
    // h_ratio_pT = make_histograms_output[1];

    // TH1D * h_signif_pT = new TH1D("h_sig_pT", "; Muon pT Cut [GeV]; Significance;", n_pT_cuts, min_cut, min_cut + n_pT_cuts*step_cut);
    // TH1D * h_ratio_pT = new TH1D("h_ratio_pT", "; Muon pT Cut [GeV]; Signal to Background Ratio;", n_pT_cuts, min_cut, min_cut + n_pT_cuts*step_cut);

    // for(int n = 0; n < n_pT_cuts; ++n){

    //     h_signal_cuts = (TH1D*) signal -> Get(h_pT_cuts[n]);
    //     h_bkgd_cuts = (TH1D*) bkgd -> Get(h_pT_cuts[n]);

    //     higgs_min = h_signal_cuts -> FindBin(100);
    //     higgs_max = h_signal_cuts -> FindBin(150);
    //     signal_data = h_signal_cuts -> Integral(higgs_min, higgs_max);
    //     bkgd_data = h_bkgd_cuts -> Integral(higgs_min, higgs_max);
    //     significance = signal_data/TMath::Sqrt(bkgd_data);
    //     ratio = signal_data/bkgd_data;

    //     std::cout << "significance = " << significance << std::endl;
    //     std::cout << "signal to background ratio = " << ratio << std::endl;

    //     bin_pt = h_signif_pT -> FindBin(pT_cut_values.at(n));

    //     bool infinity_sig = isinf(significance);
    //     bool nan_sig = isnan(significance);
    //     bool infinity_ratio = isinf(ratio);
    //     bool nan_ratio = isnan(ratio);

    //     if(!infinity_sig && !nan_sig && !infinity_ratio && !nan_ratio){
    //         h_signif_pT -> SetBinContent(bin_pt, significance);
    //         h_ratio_pT -> SetBinContent(bin_pt, ratio);
    //     } 
            
    // }

    // OutputFile->cd();
    // h_signif->Write();
    // h_ratio->Write();

    OutputFile -> Close();

}

// void Make_Histograms(TH1D * h_signif, TH1D * h_ratio, TString property, TString units, int n_cuts, double min_cut, double step_cut){
//     std::vector<TH1D*> output;
//     h_signif = new TH1D("h_signif_"+property, "; Muon" + property + "Cut [" + units + "]; Significance;", n_cuts, min_cut, min_cut + n_cuts*step_cut);
//     h_ratio = new TH1D("h_ratio_"+property, "; Muon" + property + "Cut [" + units + "]; Signal to Background Ratio;", n_cuts, min_cut, min_cut + n_cuts*step_cut);
// }

// std::vector<double> Calculate_Parameters(TH1D * h_signal, TH1D * h_bkgd, std::vector<TString> h_varycuts, int n_cuts){
    // std::vector<double> output;
    // Int_t higgs_min;
    // Int_t higgs_max;
    // double signal_data;
    // double bkgd_data;
    // double significance;
    // double ratio;
    // int bin_pt;

    // for(int n = 0; n < n_cuts; ++n){
    //     h_signal = (TH1D*) signal -> Get(h_varycuts[n]);
    //     h_bkgd = (TH1D*) bkgd -> Get(h_varycuts[n]);

    //     higgs_min = h_signal -> FindBin(100);
    //     higgs_max = h_signal -> FindBin(150);
    //     signal_data = h_signal -> Integral(higgs_min, higgs_max);
    //     bkgd_data = h_bkgd -> Integral(higgs_min, higgs_max);
    //     significance = signal_data/TMath::Sqrt(bkgd_data);
    //     ratio = signal_data/bkgd_data;

        // bin_pt = h_signif -> FindBin(pT_cut_values.at(n));

        // bool infinity_sig = isinf(significance);
        // bool nan_sig = isnan(significance);
        // bool infinity_ratio = isinf(ratio);
        // bool nan_ratio = isnan(ratio);
//     }
// }

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

        higgs_min = h_signal_cuts -> FindBin(100);
        higgs_max = h_signal_cuts -> FindBin(150);
        signal_data = h_signal_cuts -> Integral(higgs_min, higgs_max);
        bkgd_data = h_bkgd_cuts -> Integral(higgs_min, higgs_max);
        significance = signal_data/TMath::Sqrt(bkgd_data);
        ratio = signal_data/bkgd_data;

        std::cout << "significance = " << significance << std::endl;
        std::cout << "signal to background ratio = " << ratio << std::endl;

        std::vector<double> cut_values = Define_Cut_Values(n_cuts, min_cut, step_cut);

        bin_pt = h_signif -> FindBin(cut_values.at(n));

        bool infinity_sig = isinf(significance);
        bool nan_sig = isnan(significance);
        bool infinity_ratio = isinf(ratio);
        bool nan_ratio = isnan(ratio);

        if(!infinity_sig && !nan_sig && !infinity_ratio && !nan_ratio){
            std::cout << "hello" << std::endl;
            h_signif -> SetBinContent(bin_pt, significance);
            h_ratio -> SetBinContent(bin_pt, ratio);
        }
    }

    OutputFile->cd();
    h_signif->Write();
    h_ratio->Write();

}


            
    