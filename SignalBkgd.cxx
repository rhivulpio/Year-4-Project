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

TFile * signal_file;
TFile * bkgd_file;
TFile * bkgd_file_two;
TFile * OutputFile;

void SignalBkgd(){
    TString signal_file_name = "signal.root";
    TString bkgd_file_name = "bkgd.root";
    TString bkgd_file_two_name = "bkgd_two.root";

    signal_file = TFile::Open(signal_file_name); //importing signal and background files
    bkgd_file = TFile::Open(bkgd_file_name);
    bkgd_file_two = TFile::Open(bkgd_file_two_name);
    OutputFile = new TFile("SignalBkgd.root", "recreate");
    OutputFile -> cd();
    TCanvas * c2 = new TCanvas("c2"," Efficiency ", 50, 50, 1000, 750);
    gStyle -> SetOptStat(0);

    TH1D * h_signal;
    TH1D * h_bkgd;
    TH1D * h_bkgd_two;
    TLegend * legend;
    Double_t histmax;
    
    std::vector<TString> histograms; //creating a vector to store histogram names for mass reconstruction
    histograms.push_back("Mass Reconstruction/h_Higgs_reco");
    histograms.push_back("Mass Reconstruction/h_Z_reco");
    histograms.push_back("Mass Reconstruction/h_Zstar_reco");

    std::vector<TString> smeared_histograms; //creating a vector to store histogram names for smeared mass reconstruction
    smeared_histograms.push_back("Smearing/h_Higgs_reco_smeared");
    smeared_histograms.push_back("Smearing/h_Z_reco_smeared");
    smeared_histograms.push_back("Smearing/h_Zstar_reco_smeared");

    TH1D * h_signal_data;
    TH1D * h_bkgd_data;

    h_signal_data = (TH1D*) signal_file -> Get(histograms[1]); //retrieving histograms from vector that stores names
    h_bkgd_data = (TH1D*) bkgd_file -> Get(histograms[1]);

    int min = h_signal_data -> FindBin(20);
    int max = h_signal_data -> FindBin(120);
    double Z_signal = h_signal_data -> GetEntries();
    double Z_bkgd = h_bkgd_data -> GetEntries();

    //------------------------------------------------------------
    //Plotting Mass Reconstruction for Higgs, Z and Z* on One Plot
    //------------------------------------------------------------
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
        legend->SetHeader("Histogram Markers", "C");
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

    //----------------------------------------------------
    //Plotting Stacked Mass Reconstruction
    //----------------------------------------------------
    hstack(smeared_histograms[0]);

    OutputFile -> Close();
}

//--------------------------------------------------------------------
//Calculating and Plotting Significance and Signal to Background Ratio
//--------------------------------------------------------------------
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

//--------------------------------------------------------------------
//Plotting Stacked Mass Reconstruction
//--------------------------------------------------------------------
TCanvas * hstack(TString histogram) {
    THStack * hs = new THStack("hs", "Stacked Mass Reconstruction");
    
    TH1D * h2 = (TH1D*) bkgd_file -> Get(histogram); 
    h2->SetFillColor(29);
    h2->SetLineColor(1);
    hs->Add(h2);

    TH1D * h3 = (TH1D*) bkgd_file_two -> Get(histogram); 
    h3->SetFillColor(46);
    h3->SetLineColor(1);
    hs->Add(h3);
    
    TH1D * h1 = (TH1D*) signal_file -> Get(histogram);
    h1->SetFillColor(38);
    h1->SetLineColor(1);
    hs->Add(h1);

    TCanvas * cst = new TCanvas("cst", "Stacked Mass Reconstruction", 50, 50, 1000, 750);
    TH1 * last_stack = (TH1*)hs->GetStack()->Last();
    last_stack->SetFillColor(38);
    last_stack->SetLineColor(1);
    hs->Draw("HIST");
    last_stack->Draw("SAME E1");

    hs->GetXaxis()->SetTitle("m_{4l} (GeV)");
    hs->GetYaxis()->SetTitle("Number of Events");

    TLegend * legend2;
    legend2 = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend2->AddEntry(h1, "Signal");
    legend2->AddEntry(h2, "ZZ* -> 4l Background");
    legend2->AddEntry(h3, "Z -> 4l Background");
    legend2->Draw("same");

    int min_signal = h1 -> FindBin(123);
    int max_signal = h1 -> FindBin(127);
    double signal_data = h1 -> Integral(min_signal, max_signal);
    std::cout << " signal counts between 123 and 127: " << signal_data << std::endl;

    int min_bkgd = h2 -> FindBin(123);
    int max_bkgd = h2 -> FindBin(127);
    double bkgd_data = h2 -> Integral(min_bkgd, max_bkgd);
    std::cout << " bkgd counts between 123 and 127: " << bkgd_data << std::endl;

    int min_bkgd_two = h3 -> FindBin(123);
    int max_bkgd_two = h3 -> FindBin(127);
    double bkgd_data_two = h3 -> Integral(min_bkgd_two, max_bkgd_two);
    std::cout << " bkgd two counts between 123 and 127: " << bkgd_data_two << std::endl;


    hs->Write();

    return cst;
}