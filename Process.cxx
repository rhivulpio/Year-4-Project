
#include "Process.h"
#include "Utilities.cxx"
#include "GlobalConstants.cxx"

bool Debug = false;

int main(int argc, char* argv[]) {

    // Input Delphes File
 
    const TString InputFile = argv[1];
    const TString OutputFileName = argv[2];

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Running Process"  << std::endl;
    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "InputFile = " << InputFile << std::endl;
    std::cout << "OutputFileName = " << OutputFileName << std::endl;
    std::cout << "-------------------------------------------------------------"  << std::endl;

    ExRootTreeReader * reader = NULL;
    reader = InitReader(InputFile);

    //------------------------------------
    // Declare the output
    //------------------------------------

    OutputFile = new TFile(OutputFileName,"recreate");

    OutputFile->cd();
    OutputFile->mkdir("Transverse Momentum");
    OutputFile->mkdir("Pseudorapidity");
    OutputFile->mkdir("Transverse Energy");
    OutputFile->mkdir("Higgs Plots");
    OutputFile->mkdir("4mu Events");
    OutputFile->mkdir("4mu Events Seen by Detector");
    OutputFile->mkdir("Acceptance Plots");
    OutputFile->mkdir("True Jets");
    OutputFile->mkdir("Missing Energy");
    OutputFile->mkdir("Scatter Plots");
    OutputFile->mkdir("Kinematic Reconstruction");
    OutputFile->mkdir("Mass Reconstruction");
    OutputFile->mkdir("Event Weight");
    OutputFile->mkdir("Vary Cuts/pT Cuts");
    OutputFile->mkdir("Vary Cuts/Zstar Cuts");
    OutputFile->mkdir("Vary Cuts/Z Cuts");
    OutputFile->mkdir("Smearing");
   

    h_EventCount = new TH1D("h_EventCount",";""; Number of Events", 9, 0, 9);
    // gStyle -> SetOptStat(0);
    h_EventCount->SetStats(0);
    TAxis * xAxis = h_EventCount->GetXaxis();
    xAxis->SetBinLabel(1, "Total Events");
    xAxis->SetBinLabel(2, "4mu Events");
    xAxis->SetBinLabel(3, "4mu Events Seen");
    xAxis->SetBinLabel(4, "4e Events");
    xAxis->SetBinLabel(5, "4e Events Seen");
    xAxis->SetBinLabel(6, "2e2mu Events");
    xAxis->SetBinLabel(7, "2e2mu Events Seen");
    xAxis->SetBinLabel(8, "4l Events");
    xAxis->SetBinLabel(9, "4l Events Seen");
    // h_EventCount->Draw("");
    // PaintBin(h_EventCount, 3, kRed);
    // PaintBin(h_EventCount, 9, kBlue);

    h_WeightCount = new TH1D("h_WeightCount","", 1, 0, 1);

    h_ZZ_Mass = new TH1D("h_ZZ_Mass","; ZZ Mass [GeV]; Events / 2 GeV", 125, 0.0, 250.0);

    //transverse momentum histograms
    h_Lepton_Pt = new TH1D("h_Lepton_Pt","; Charged Lepton p_{T} [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h_Z_Pt = new TH1D("h_Z_Pt","; Z p_{T} [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h_Jet_Pt = new TH1D("h_Jet_Pt","; Jet p_{T} [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h_mu_pT = new TH1D("h_mu_pT", "; Muon Transverse Momentum [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h_Higgs_pT = new TH1D("h_Higgs_pT","; Higgs Transverse Momentum ; Events", 100, 0, 170);

    //pseudorapidity histograms
    h_Jet_eta = new TH1D("h_Jet_eta","; Jet Eta ; Events", 50, -6.0, 1.0);
    h_Z_eta = new TH1D("h_Z_eta","; Z Boson Eta ; Events", 50, -7.0, 1.0);
    h_mu_eta = new TH1D("h_mu_eta","; Muon Eta ; Events", 50, -7.0, 5.0);
    h_nu_eta = new TH1D("h_nu_eta","; Neutrino Eta ; Events", 50, -6.0, 4.0);
    h_Higgs_eta = new TH1D("h_Higgs_eta","; Higgs Eta ; Events", 100, -8.0, 0); 

    //transverse energy histograms
    h_Jet_Et = new TH1D("h_Jet_Et","; Jet Transverse Energy [GeV]; Events", 100, 0, 160);
    h_Z_Et = new TH1D("h_Z_Et","; Z Transverse Energy [GeV]; Events", 100, 0, 180);
    h_mu_Et = new TH1D("h_mu_Et", "; Muon Transverse Energy [GeV] ; Events", 100, 0, 140);
    h_nu_Et = new TH1D("h_nu_Et","; Neutrino Transverse Energy [GeV] ; Events", 100, 0, 200);
    h_Higgs_Et = new TH1D("h_Higgs_Et","; Higgs Transverse Energy [GeV]; Events", 100, 100, 240);
    
    //4mu event histograms
    h4mu_mu_pT = new TH1D("h4mu_mu_pT", "; Muon Transverse Momentum [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h4mu_mu_eta = new TH1D("h4mu_mu_eta","; Muon Eta ; Events", 50, -7.0, 5.0);
    h4mu_nu_eta = new TH1D("h4mu_nu_eta","; Neutrino Eta ; Events", 50, -6.0, 4.0);
    h4mu_mu_Et = new TH1D("h4mu_mu_Et", "; Muon Transverse Energy [GeV] ; Events", 100, 0, 140);
    h4mu_nu_pT = new TH1D("h4mu_nu_pT","; Neutrino Transverse Momentum [GeV] ; Events", 100, 0, 200);

    //4mu event seen by detector histograms
    h4mu_mu_pT_seen = new TH1D("h4mu_mu_pT_seen", "; Muon Transverse Momentum [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h4mu_mu_eta_seen = new TH1D("h4mu_mu_eta_seen","; Muon Eta ; Events", 50, -7.0, 5.0);
    h4mu_nu_eta_seen = new TH1D("h4mu_nu_eta_seen","; Neutrino Eta ; Events", 50, -5.0, 5.0);
    h4mu_mu_Et_seen = new TH1D("h4mu_mu_Et_seen", "; Muon Transverse Energy [GeV] ; Events", 100, 0, 140);
    h4mu_nu_pT_seen = new TH1D("h4mu_nu_pT_seen","; Neutrino Transverse Momentum [GeV] ; Events", 100, 0, 200);
    h4mu_Z_pT_seen = new TH1D("h4mu_Z_pT_seen","; Z Boson Transverse Momentum [GeV] ; Events", 100, 0, 200);
    h4mu_Z_eta_seen = new TH1D("h4mu_Z_eta_seen","; Z Boson Pseudorapidity ; Events", 100, -7, 2);
    h4mu_Higgs_pT_seen = new TH1D("h4mu_Higgs_pT_seen","; Higgs Boson Transverse Momentum [GeV] ; Events", 100, 0, 300);
    h4mu_Higgs_eta_seen = new TH1D("h4mu_Higgs_eta_seen","; Higgs Boson Pseudorapidity ; Events", 100, -8, 0);
    h4mu_jet_eta_seen = new TH1D("h4mu_jet_eta_seen","; Jet Pseudorapidity ; Events", 100, -7, 1);
    h4mu_jet_pT_seen = new TH1D("h4mu_jet_pT_seen","; Jet Transverse Momentum [GeV] ; Events", 100, 0, 200);

    h4mu_mu_pT = new TH1D("h4mu_mu_pT", "; Muon Transverse Momentum [GeV]; Events / 5 GeV", 200, 0.0, 300.0);
    h4mu_mu_eta = new TH1D("h4mu_mu_eta","; Muon Eta ; Events", 50, -7.0, 5.0);
    h4mu_nu_eta = new TH1D("h4mu_nu_eta","; Neutrino Eta ; Events", 50, -6.0, 4.0);
    h4mu_nu_pT = new TH1D("h4mu_nu_pT","; Neutrino Transverse Momentum [GeV] ; Events", 100, 0, 200);

    //Acceptance plots
        //for all events
    e_eta = new TEfficiency("e_eta", "Plot of Acceptance against Pseudorapidity of Muons; Muon Pseudorapidity; Acceptance", 100, -10, 10);
    e_Et = new TEfficiency("e_Et", "Plot of Acceptance against Transverse Energy of Muons; Muon Transverse Energy; Acceptance", 100, -100, 100);
    e_pT = new TEfficiency("e_pT", "Plot of Acceptance against Transverse Momentum of Muons; Muon Transverse Momentum; Acceptance", 100, -100, 100);
    e_H_eta = new TEfficiency("e_H_eta", "Plot of Acceptance against Pseudorapidity of Higgs (cuts); Higgs Pseudorapidity; Acceptance", 100, -10, 10);
    e_H_Et = new TEfficiency("e_H_Et", "Plot of Acceptance against Transverse Energy of Higgs (cuts); Higgs Transverse Energy; Acceptance", 100, -100, 100);
    e_H_pT = new TEfficiency("e_H_pT", "Plot of Acceptance against Transverse Momentum of Higgs (cuts); Higgs Transverse Momentum; Acceptance", 100, -100, 100);
        //for 4mu events only
    e4mu_eta = new TEfficiency("e4mu_eta", "Plot of Acceptance against Pseudorapidity of Muons; Muon Pseudorapidity; Acceptance", 100, -10, 10);
    e4mu_Et = new TEfficiency("e4mu_Et", "Plot of Acceptance against Transverse Energy of Muons; Muon Transverse Energy; Acceptance", 100, -100, 100);
    e4mu_pT = new TEfficiency("e4mu_pT", "Plot of Acceptance against Transverse Momentum of Muons; Muon Transverse Momentum; Acceptance", 100, -100, 100);
    e4mu_H_eta = new TEfficiency("e4mu_H_eta", "Plot of Acceptance against Pseudorapidity of Higgs; Higgs Pseudorapidity; Acceptance", 100, -10, 10);
    e4mu_H_Et = new TEfficiency("e4mu_H_Et", "Plot of Acceptance against Transverse Energy of Higgs; Higgs Transverse Energy; Acceptance", 100, -100, 100);
    e4mu_H_pT = new TEfficiency("e4mu_H_pT", "Plot of Acceptance against Transverse Momentum of Higgs; Higgs Transverse Momentum; Acceptance", 100, -100, 100);

    //Plotting different variables against each other - need to look at this
    e_mu_pT_eta = new TEfficiency("e_mu_pT_eta", "; Muon Pseudorapidity; Muon Transverse Momentum [GeV]", 100, -10, 10);
    h_mu_pT_eta = new TH2D("h_mu_pT_eta", "Muon Transverse Momentum against Muon Pseudorapidity; Muon Pseudorapidity; Muon Transverse Momentum [GeV]", 30, -8, 4, 30, 0, 180);
    h_4mu_pT_eta = new TH2D("h_4mu_pT_eta", "Muon Transverse Momentum against Muon Pseudorapidity for 4mu events; Muon Pseudorapidity; Muon Transverse Momentum [GeV]", 30, -8, 4, 30, 0, 180);

    //true jet histograms
    h_trueJet_Et = new TH1D("h_trueJet_Et","; True Jet Transverse Energy [GeV]; Events", 200, 0, 180);
    h_trueJet_eta = new TH1D("h_trueJet_eta","; True Jet Pseudorapidity; Events", 200, -6, 0);
    h_trueJet_Pt = new TH1D("h_trueJet_Pt","; True Jet Transverse Momentum [GeV]; Events", 200, 0, 200);

    //2D plot of missing energy against neutrino energy - should be linear
    h_ME_nu_pT = new TH2D("h_ME_nu_pT","Missing Transverse Momentum against Neutrino Transverse Momentum; Neutrino Transverse Momentum [GeV]; Missing Transverse Momentum [GeV]", 100, 0, 200, 100, 0, 200);
    h_ME_nu_eta = new TH2D("h_ME_nu_eta","Missing Eta against Neutrino Eta; Neutrino Eta; Missing Eta", 100, -4, 5, 100, -7, 5);
    h_ME_nu_phi = new TH2D("h_ME_nu_phi","Missing Phi against Neutrino Phi; Neutrino Phi; Missing Phi", 100, -4, 4, 100, -4, 4);

    //kinematic reconstruction plots
        //electron reconstruction method
    x_Qsquared_electron = new TH2D("x_Qsquared_electron", "Q^{2} against x (Electron Reconstruction Method); x; Q^{2}", 100, 0, 0.12, 100, 0, 20000);
    h_logx_electron = new TH1D("h_logx_electron", "; log_{10}x; Events", 100, -6, 1);
    h_logy_electron = new TH1D("h_logy_electron", "; log_{10}y; Events", 100, -2, 0.5);
    h_logQsquared_electron = new TH1D("h_logQsquared_electron", "; log_{10}Q^{2}; Events", 100, -1, 6);
        //hadron reconstruction method
    x_Qsquared_hadron = new TH2D("x_Qsquared_hadron", "Q^{2} against x (Hadron Reconstruction Method); x; Q^{2}", 100, 0, 0.12, 100, 0, 20000);
    h_logx_hadron = new TH1D("h_logx_hadron", "; log_{10}x; Events", 100, -6, 0);
    h_logy_hadron = new TH1D("h_logy_hadron", "; log_{10}y; Events", 100, -2, 0.5);
    h_logQsquared_hadron = new TH1D("h_logQsquared_hadron", "; log_{10}Q^{2}; Events", 100, -1, 6);
    h_logx_hadron->SetStats(0);
    h_logy_hadron->SetStats(0);
    h_logQsquared_hadron->SetStats(0);
        //comparing electron and hadron reconstruction methods
    log_Qsquared_plot = new TH2D("log_Qsquared_plot", "log_{10}(Q^{2}) (Electron) against log_{10}(Q^{2}) (Hadron); log_{10}(Q^{2}_{h}); log_{10}(Q^{2}_{e})", 100, 0, 7, 100, 0, 7);
    log_x_plot = new TH2D("log_x_plot", "log_{10}(x) (Electron) against log_{10}(x) (Hadron); log(x_{h}); log_{10}(x_{e})", 100, -5, 0, 100, -5, 0);
    log_y_plot = new TH2D("log_y_plot", "log_{10}(y) (Electron) against log_{10}(y) (Hadron); log(y_{h}); log_{10}(y_{e})", 100, -2, 1, 100, -2, 1);

    //mass reconstruction
    h_Higgs_reco = new TH1D("h_Higgs_reco", "; m_{4l} [GeV]; Events ", 150, 80.0, 200.0);
    h_ZZ_mass_reco = new TH1D("h_ZZ_mass_reco", "; Reconstructed ZZ* Mass [GeV]; Events ", 150, 0.0, 200.0);
    h_Z_reco = new TH1D("h_Z_reco", "; m_{ll} (Leading Lepton Pair) [GeV]; Events ", 150, 0.0, 200.0);
    h_Zstar_reco = new TH1D("h_Zstar_reco", "; m_{ll} (Subleading Lepton Pair) [GeV]; Events ", 150, 0.0, 200.0);

    //event weight plot
    h_eventweight = new TH1D("h_eventweight", "; Event Weight ; Events ", 100, 0.000100, 0.000120);

    //smeared mass reconstruction
    h_Higgs_reco_smeared = new TH1D("h_Higgs_reco_smeared", "; m_{4l} [GeV]; Events ", 150, 80.0, 200.0);
    //h_ZZ_mass_reco_smeared = new TH1D("h_ZZ_mass_reco_smeared", "; Reconstructed ZZ* Mass [GeV]; Events ", 150, 0.0, 130.0);
    h_Z_reco_smeared = new TH1D("h_Z_reco_smeared", "; m_{ll} (Leading Lepton Pair) [GeV]; Events ", 150, 0.0, 200.0);
    h_Zstar_reco_smeared = new TH1D("h_Zstar_reco_smeared", "; m_{ll} (Subleading Lepton Pair) [GeV]; Events ", 150, 0.0, 200.0);

    // Run the selection
    Process(reader);

    std::cout << "Total Number of Events: " << h_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Number of 4mu Events: " << h_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Number of 4mu Events Seen by Detector: " << h_EventCount->GetBinContent(3) << std::endl; 
    std::cout << "Number of 4e Events: " << h_EventCount->GetBinContent(4) << std::endl;
    std::cout << "Number of 4e Events Seen by Detector: " << h_EventCount->GetBinContent(5) << std::endl;
    std::cout << "Number of 2e2mu Events: " << h_EventCount->GetBinContent(6) << std::endl;
    std::cout << "Number of 2e2mu Events Seen by Detector: " << h_EventCount->GetBinContent(7) << std::endl;
    std::cout << "Number of 4l Events: " << h_EventCount->GetBinContent(8) << std::endl;
    std::cout << "Number of 4l Events Seen by Detector: " << h_EventCount->GetBinContent(9) << std::endl;

    std::cout << "Write to file..." << std::endl;

    OutputFile->cd(); 

    h_EventCount->Write();

    h_WeightCount->Write();

    OutputFile->cd("Transverse Momentum");

    h_Z_Pt->Write();
    h_Lepton_Pt->Write();
    h_Jet_Pt->Write();
    h_ZZ_Mass->Write();
    h_mu_pT->Write();

    OutputFile->cd("Pseudorapidity");

    h_mu_eta->Write();
    h_Jet_eta->Write();
    h_Z_eta->Write();
    h_nu_eta->Write();

    OutputFile->cd("Transverse Energy");

    h_mu_Et->Write();
    h_Jet_Et->Write();
    h_Z_Et->Write();
    h_nu_Et->Write();

    OutputFile->cd("Higgs Plots");

    h_Higgs_pT->Write();
    h_Higgs_eta->Write();
    h_Higgs_Et->Write();

    OutputFile->cd("4mu Events");

    h4mu_mu_pT->Write();
    h4mu_mu_Et->Write();
    h4mu_nu_pT->Write();
    h4mu_mu_eta->Write();
    h4mu_nu_eta->Write(); 

    OutputFile->cd("4mu Events Seen by Detector");
    h4mu_mu_pT_seen->Write();
    h4mu_mu_Et_seen->Write();
    h4mu_nu_pT_seen->Write();
    h4mu_mu_eta_seen->Write();
    h4mu_nu_eta_seen->Write(); 
    h4mu_Z_pT_seen->Write();
    h4mu_Z_eta_seen->Write();
    h4mu_Higgs_pT_seen->Write();
    h4mu_Higgs_eta_seen->Write();
    h4mu_jet_eta_seen->Write();
    h4mu_jet_pT_seen->Write();

    OutputFile->cd("Acceptance Plots");

    e_eta->Write();
    e_Et->Write();
    e_pT->Write();
    e4mu_eta->Write();
    e4mu_Et->Write();
    e4mu_pT->Write();

    e_H_eta->Write();
    e_H_Et->Write();
    e_H_pT->Write();
    e4mu_H_eta->Write();
    e4mu_H_Et->Write();
    e4mu_H_pT->Write();

    e_mu_pT_eta->Write();

    OutputFile->cd("True Jets");

    h_trueJet_Et->Write();
    h_trueJet_eta->Write();
    h_trueJet_Pt->Write();

    OutputFile->cd("Missing Energy");

    h_ME_nu_pT->Write();
    h_ME_nu_eta->Write();
    h_ME_nu_phi ->Write();

    OutputFile->cd("Scatter Plots");

    h_mu_pT_eta->Write();
    h_4mu_pT_eta->Write();

    OutputFile->cd("Kinematic Reconstruction");

    x_Qsquared_electron->Write();
    x_Qsquared_hadron->Write();
    log_Qsquared_plot->Write();
    log_x_plot->Write();
    log_y_plot->Write();
    h_logx_electron->Write();
    h_logy_electron->Write();
    h_logQsquared_electron->Write();
    h_logx_hadron->Write();
    h_logy_hadron->Write();
    h_logQsquared_hadron->Write();

    OutputFile->cd("Mass Reconstruction");
    h_Higgs_reco->Write();
    h_ZZ_mass_reco->Write();
    h_Z_reco->Write();
    h_Zstar_reco->Write();

    TCanvas * c1 = new TCanvas("c1", "Mass Reconstruction", 50, 50, 1500, 1200);
    TLegend * legend;

    h_Higgs_reco->SetTitle("Reconstruction of Higgs, Z and Z* Masses");
    h_Higgs_reco->GetXaxis()->SetTitle("Mass (GeV)");
    h_Higgs_reco->GetYaxis()->SetTitle("Number of Events");
    h_Higgs_reco->GetXaxis()->SetRangeUser(0, 200);
    h_Higgs_reco->SetLineColor(kRed);
    h_Higgs_reco->SetStats(kFALSE);
    h_Higgs_reco->Draw("hist E2");

    // h_Z_reco->SetTitle("Reconstruction of Higgs, Z and Z* Masses");
    // h_Z_reco->GetXaxis()->SetTitle("Mass (GeV)");
    // h_Z_reco->GetYaxis()->SetTitle("Number of Events");
    // h_Z_reco->GetXaxis()->SetRangeUser(0, 130);
    // h_Z_reco->SetLineColor(kRed);
    // h_Z_reco->SetStats(kFALSE);
    // h_Z_reco->Draw("hist E2");

    h_Z_reco->SetLineColor(kGreen);
    h_Z_reco->Draw("hist same E2");

    h_Zstar_reco->SetLineColor(kBlue);
    h_Zstar_reco->Draw("hist same E2");

    legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->SetHeader("Particle", "C");
    legend->AddEntry(h_Higgs_reco, "Higgs Boson");
    legend->AddEntry(h_Z_reco, "Z Boson");
    legend->AddEntry(h_Zstar_reco, "Z* Boson");
    legend->Draw("same");

    c1->Write("Mass Reconstruction");

    OutputFile->cd("Event Weight");

    h_eventweight->Write();

    OutputFile->cd("Vary Cuts/pT Cuts");
    Write_Histogram(h_pT_cuts);

    OutputFile->cd("Vary Cuts/Zstar Cuts");
    Write_Histogram(h_Zstar_cuts);

    OutputFile->cd("Vary Cuts/Z Cuts");
    Write_Histogram(h_Z_cuts);

    OutputFile->cd("Smearing");
    h_Higgs_reco_smeared->Write();
    h_Z_reco_smeared->Write();
    h_Zstar_reco_smeared->Write();

    OutputFile->Close();

    std::cout << "Tidy..." << std::endl;

    delete reader;

    std::cout << "Done!" << std::endl;

    return 0;

}

ExRootTreeReader * InitReader(const TString FilePath) {

    std::cout << "InitReader" << std::endl;

    TFile * f = TFile::Open(FilePath);

    TChain * Chain = new TChain("Delphes","");

    Chain->Add(FilePath);

    // Create object of class ExRootTreeReader
    ExRootTreeReader * r = new ExRootTreeReader(Chain);

    return r;
}

void Process(ExRootTreeReader * treeReader) {

    // Get pointers to branches used in this analysis
    bEvent = treeReader->UseBranch("Event");
    bJet = treeReader->UseBranch("GenJet");
    bTruthLepton = treeReader->UseBranch("TruthLeptonParticles");
    bTruthWZ = treeReader->UseBranch("TruthWZHParticles");

    Long64_t numberOfEntries = treeReader->GetEntries();

    if (Debug) numberOfEntries = 1000;

    int nSelected = 0;

    std::cout << "-------------------------------------------------------------"  << std::endl;
    std::cout << "Input: " << numberOfEntries << " events to process" << std::endl;
    
    //------------------------------------------------------------------
    //calculating scale factors for signal and background files
    //------------------------------------------------------------------
    double sigma_signal = 1.34e-17; //in picobarns
    double sigma_bkgd = 8.9e-18; //in picobarns
    double sigma_bkgd_two = 2.41e-18; //in picobarns
    double no_events = 100000;
    double no_events_two = 10000;
    double luminosity_signal = no_events/sigma_signal;
    double luminosity_bkgd = no_events/sigma_bkgd;
    double luminosity_bkgd_two = no_events_two/sigma_bkgd_two;
    double luminosity_LHeC = 1000e15; //1000 inverse femtobarns
    double sf_signal = luminosity_LHeC/luminosity_signal; //scale factor
    double sf_bkgd = luminosity_LHeC/luminosity_bkgd;
    double sf_bkgd_two = luminosity_LHeC/luminosity_bkgd_two;

    std::cout << "sf signal = " << sf_signal << " sf bkgd = " << sf_bkgd << " sf bkgd 2 = " << sf_bkgd_two << std::endl;

    std::vector<bool> pT_cuts; //defines an array which will contain a flag to indicate whether each event passes each minimum pT cut

    h_pT_cuts = Define_Histograms("h_pT", n_pT_cuts);
    h_Zstar_cuts = Define_Histograms("h_Zstar", n_Zstar_cuts);
    h_Z_cuts = Define_Histograms("h_Z", n_Z_cuts);

    std::vector<double> pT_cut_values = Define_Cut_Values(n_pT_cuts, min_pT_cut, step_pT_cut);
    std::vector<double> Zstar_cut_values = Define_Cut_Values(n_Zstar_cuts, min_Zstar_cut, step_Zstar_cut);
    std::vector<double> Z_cut_values = Define_Cut_Values(n_Z_cuts, min_Z_cut, step_Z_cut);

    //------------------------------------------------------------------
    // Loop Over All Events
    //------------------------------------------------------------------
    for(Int_t entry = 0; entry < numberOfEntries; ++entry){

        //Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0);

        //const float Event_Weight = signal ? sf_signal : sf_bkgd;
        const double event_weight = sf_signal; //manually change this

        h_EventCount->Fill(0.5, event_weight);
        h_WeightCount->Fill(0.5, event_weight);
        h_eventweight->Fill(event_weight);

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //This is to calculate the missing energy, set it equal to zero at the start of each event and the vector 
        //sum will accumulate throughout.
        TLorentzVector missing_energy_vector;
        missing_energy_vector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);

        //------------------------------------------------------------------
        // Jet Loop
        //------------------------------------------------------------------
        bool jet_cut = true; //flag to check if the jet passes the jet cut
        double eta_min_j = -5.0; //eta range of the hadronic calorimeter
        double eta_max_j = 5.5;
        double pT_min = 15; //in GeV

        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet * jet = (Jet*) bJet->At(i);

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << 
            " mass = " << jet->Mass << " flavour = " << jet->Flavor << std::endl;

            TLorentzVector vec_jet;
            vec_jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
            
            //defining the selection criteria for the jets
            if(
                eta_min_j>jet->Eta && 
                jet->Eta>eta_max_j && 
                jet->PT<pT_min){
                jet_cut = false;
            }

            if(jet_cut){
                h_Jet_Pt->Fill(vec_jet.Pt(), event_weight);
                h_Jet_eta->Fill(vec_jet.Eta(), event_weight);
                h_Jet_Et->Fill(TMath::Sqrt(vec_jet.Pt() * vec_jet.Pt() + vec_jet.M() * vec_jet.M()), event_weight);
            }

            TLorentzVector vec_lep;
            bool true_jet = true; //flag to check if the jet is a true jet or not

            //lepton loop inside jet loop to calculate DeltaR for each jet/lepton combination
            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) { 

                GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

                vec_lep.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);

                double delta_R = TMath::Power(vec_jet.Phi() - vec_lep.Phi(), 2) + TMath::Power(vec_jet.Eta() - vec_lep.Eta(), 2);
                double jet_radius = 0.4; //defined in the DELPHES configuration

                //a true jet is defined as one where the DeltaR between the jet and the lepton is greater than the jet radius
                if(delta_R < jet_radius) true_jet = false;
            }

            if(Debug) std::cout << "  Is it a true jet? " << (true_jet ? "yes" : "no") << std::endl;
        
            //------------------------------------------------------------------
            // Jet Loop for True Jets That Pass the Jet Cuts
            //------------------------------------------------------------------
            if(true_jet && jet_cut){
                //adding true jet four vectors to the missing energy four vector
                missing_energy_vector = missing_energy_vector + vec_jet; 
                if(Debug) std::cout << "  Missing Energy Vector (Jets) " << missing_energy_vector.Pt() << std::endl;

                h_trueJet_Pt->Fill(vec_jet.Pt(), event_weight);
                h_trueJet_eta->Fill(vec_jet.Eta(), event_weight);
                h_trueJet_Et->Fill(TMath::Sqrt(vec_jet.Pt() * vec_jet.Pt() + vec_jet.M() * vec_jet.M()), event_weight);
                
            }
            
        }
        // Jet Loop End

        bool event4mu = false; //flag to check if the event satisfies the 4mu subchannel
        bool event4mu_seen;
        bool event4e = false;
        bool event4e_seen;
        bool event2e2mu = false;
        bool event2e2mu_seen;
        bool event2pos2neg = false; //flag to check if there are 2 positively charged leptons and 2 negatively charged leptons in the event

        bool muon_cut = true; //flag to check if muons pass cuts
        double eta_min_mu = -4.6; //eta range of the inner tracker
        double eta_max_mu = 5.3;
        //momentum cut has been previously defined
        bool electron_cut = true;
        double eta_min_e = -5.0;
        double eta_max_e = 5.0;

        std::vector<bool> pT_cut_flags = Initialise_Flags(n_pT_cuts);
        std::vector<bool> Zstar_cut_flags = Initialise_Flags(n_Zstar_cuts);
        std::vector<bool> Z_cut_flags = Initialise_Flags(n_Z_cuts);

        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------
        std::vector<GenParticle*> all_muons;
        std::vector<GenParticle*> all_electrons;
        std::vector<GenParticle*> all_leptons;
        std::vector<TLorentzVector> particles;
        std::vector<TLorentzVector> antiparticles;
        TLorentzVector temp_vector;
        
        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

            if(abs(lep->PID) == 11){ 
                all_electrons.push_back(lep);
            }

            if(abs(lep->PID) == 13){ 
                all_muons.push_back(lep);
            }

            all_leptons.push_back(lep); 
        }

        for(int i = 0; i < all_leptons.size(); ++i){
            temp_vector.SetPtEtaPhiM(all_leptons[i]->PT, all_leptons[i]->Eta, all_leptons[i]->Phi, all_leptons[i]->Mass);
            if(all_leptons[i]->Charge == 1){
                particles.push_back(temp_vector);
            }
            else if(all_leptons[i]->Charge == -1){
                antiparticles.push_back(temp_vector);
            }
        } 

        if(Debug){
            std::cout << " all_electrons.size = " << all_electrons.size() << std::endl;
            std::cout << " all_muons.size = " << all_muons.size() << std::endl;
            std::cout << " all_leptons.size = " << all_leptons.size() << std::endl;
            std::cout << " particles.size = " << particles.size() << std::endl;
            std::cout << " antiparticles.size = " << antiparticles.size() << std::endl;
        }

        if(particles.size() == 2 && antiparticles.size() == 2){
            event2pos2neg = true;
        }

        if(all_leptons.size() == 5 && event2pos2neg){

            if(all_electrons.size() == 4 && all_muons.size() == 0){
                event4e = true;
            }

            if(all_muons.size() == 4 && all_electrons.size() == 0){
                event4mu = true;
            }

            if(all_electrons.size() == 2 && all_muons.size() == 2){
                event2e2mu = true;
            }
        }

        if(event4mu) event4mu_seen = true;
        if(event4e) event4e_seen = true;
        if(event2e2mu) event2e2mu_seen = true;

        if(all_leptons.size() != 5){
            event4mu_seen = false;
            event4e_seen = false;
            event2e2mu_seen = false;
        }

        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

            if(Debug){
                std::cout << "Lepton " << i << " PID = " << lep->PID << " Charge = " << lep->Charge << " pT = " << lep->PT << " eta = " << lep->Eta
                << " phi = " << lep->Phi << " mass = " << lep->Mass << " energy = " << lep->E << " momentum = " << 
                lep->P << std::endl;
            }

            //------------------------------------------------------------------
            // Defining Selection Critera for Muons and Electrons
            //------------------------------------------------------------------
            if(abs(lep->PID) == 13){
                if(
                    lep->Eta > eta_min_mu &&
                    lep->Eta < eta_max_mu && 
                    lep->PT > pT_min) {
                        muon_cut = true;
                    } else{
                        muon_cut = false;
                    }
                
                if(Debug) std::cout << "  Does it pass the muon cut? " << (muon_cut ? "yes" : "no") << std::endl;

                if(!muon_cut){
                    event4mu_seen = false;
                    event2e2mu_seen = false;
                }

                pT_cut_flags = Check_Cuts(pT_cut_values, lep->PT, pT_cut_flags);
            }

            if(event4e || event2e2mu){
                event4mu_seen = false;
            }

            if(abs(lep->PID) == 11){
                if(
                    lep->Eta > eta_min_e &&
                    lep->Eta < eta_max_e && 
                    lep->PT > pT_min) {
                        electron_cut = true;
                    } else{
                        electron_cut = false;
                    }
                
                if(Debug) std::cout << "  Does it pass the electron cut? " << (electron_cut ? "yes" : "no") << std::endl;

                if(!electron_cut){
                    event4e_seen = false;
                    event2e2mu_seen = false;
                }
            }

            if(event4mu || event2e2mu){
                event4e_seen = false;
            }

            if(event4e || event4mu){
                event2e2mu_seen = false;
            }
        
            TLorentzVector vec_lepton;
            vec_lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
            
            // Look for electrons or muons
            if(abs(lep->PID) == 11 || abs(lep->PID) == 13 ){
                    h_Lepton_Pt->Fill(vec_lepton.Pt(), event_weight);  
            }

            //plotting histograms for muons that pass the cuts
            if(abs(lep->PID) == 13 && muon_cut){
                h_mu_eta -> Fill(vec_lepton.Eta(), event_weight);
                h_mu_Et -> Fill (TMath::Sqrt(vec_lepton.Pt() * vec_lepton.Pt() + vec_lepton.M() * vec_lepton.M()), event_weight);
                h_mu_pT -> Fill(vec_lepton.Pt(), event_weight);
            }

            //plotting acceptance against different variables for muons in all events
            if(abs(lep->PID) == 13){
                e_eta->FillWeighted(muon_cut, event_weight, vec_lepton.Eta());
                e_Et->FillWeighted(muon_cut, event_weight, TMath::Sqrt(vec_lepton.Pt() * vec_lepton.Pt() + vec_lepton.M() * vec_lepton.M()));
                e_pT->FillWeighted(muon_cut, event_weight, vec_lepton.Pt());
                e_mu_pT_eta->FillWeighted(vec_lepton.Pt(), event_weight, vec_lepton.Eta());
            }

            //looking for neutrinos
            if(abs(lep->PID) == 12){
                h_nu_Et -> Fill(TMath::Sqrt(vec_lepton.Pt() * vec_lepton.Pt() + vec_lepton.M() * vec_lepton.M()), event_weight);
                h_nu_eta -> Fill(vec_lepton.Eta(), event_weight);
            }
        } // Lepton Loop End

        if(event4mu){
            h_EventCount -> Fill(1.5, event_weight);
        }

        if(event4mu_seen){
            h_EventCount -> Fill(2.5, event_weight);
        }

        if(event4e){
            h_EventCount -> Fill(3.5, event_weight);
        }

        if(event4e_seen){
            h_EventCount -> Fill(4.5, event_weight);
        }

        if(event2e2mu){
            h_EventCount -> Fill(5.5, event_weight);
        }

        if(event2e2mu_seen){
            h_EventCount -> Fill(6.5, event_weight);
        }

        if(event4mu || event4e || event2e2mu){
            h_EventCount -> Fill(7.5, event_weight);
        }

        if(event4mu_seen || event4e_seen || event2e2mu_seen){
            h_EventCount -> Fill(8.5, event_weight);
        }

        if(Debug){
            std::cout << "  Is it a 4mu event? " << (event4mu ? "yes" : "no") << std::endl;
            std::cout << "  Is it a 4mu event seen by the detector? " << (event4mu_seen ? "yes" : "no") << std::endl;
            std::cout << "  Is it a 4e event? " << (event4e ? "yes" : "no") << std::endl;
            std::cout << "  Is it a 4e event seen by the detector? " << (event4e_seen ? "yes" : "no") << std::endl;
            std::cout << "  Is it a 2e2mu event? " << (event2e2mu ? "yes" : "no") << std::endl;
            std::cout << "  Is it a 2e2mu event seen by the detector? " << (event2e2mu_seen ? "yes" : "no") << std::endl;
        }
        
        //------------------------------------------------------------------
        // Loop For All 4l Events
        //------------------------------------------------------------------
        if(event4mu || event4e || event2e2mu){

            // if(Debug){
            //     std::cout << "  Is it a 4mu event? " << (event4mu ? "yes" : "no") << std::endl;
            //     std::cout << "  Is it a 4mu event seen by the detector? " << (event4mu_seen ? "yes" : "no") << std::endl;
            //     std::cout << "  Is it a 4e event? " << (event4e ? "yes" : "no") << std::endl;
            //     std::cout << "  Is it a 4e event seen by the detector? " << (event4e_seen ? "yes" : "no") << std::endl;
            //     std::cout << "  Is it a 2e2mu event? " << (event2e2mu ? "yes" : "no") << std::endl;
            //     std::cout << "  Is it a 2e2mu event seen by the detector? " << (event2e2mu_seen ? "yes" : "no") << std::endl;
            // }
            
            TLorentzVector vec_neutrino;
            std::vector<GenParticle*> all_muons_seen; //makes an array which will contain all muons that are seen by the detector
            std::vector<GenParticle*> all_electrons_seen;

            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

                GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

                TLorentzVector vec_lepton;
                vec_lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);

                //adding all of the muon and electron four vectors to the missing energy four vector
                if(lep->Status == 1 && abs(lep->PID) == 13){
                    missing_energy_vector = missing_energy_vector + vec_lepton;
                }

                if(lep->Status == 1 && abs(lep->PID) == 11){
                    missing_energy_vector = missing_energy_vector + vec_lepton;
                }

                vec_neutrino = vec_lepton;

                if(event4mu){
                    //plotting histograms for muons in 4mu events
                    if(abs(lep->PID) == 13){
                        h4mu_mu_eta -> Fill(vec_lepton.Eta(), event_weight);
                        h4mu_mu_Et -> Fill (TMath::Sqrt(vec_lepton.Pt() * vec_lepton.Pt() + vec_lepton.M() * vec_lepton.M()), event_weight);
                        h4mu_mu_pT-> Fill(vec_lepton.Pt(), event_weight);
                        h_mu_pT_eta -> Fill(vec_lepton.Eta(), vec_lepton.Pt());
                    }

                    //plotting histograms for neutrinos in 4mu events
                    if(abs(lep->PID) == 12){
                        h4mu_nu_pT -> Fill(vec_lepton.Pt(), event_weight);
                        h4mu_nu_eta -> Fill(vec_lepton.Eta(), event_weight); 
                    }
                }
                
                //plots for 4mu events that are seen by the detector
                if(event4mu_seen){
                    if(abs(lep->PID) == 13){
                        all_muons_seen.push_back(lep);
                        h4mu_mu_pT_seen->Fill(vec_lepton.Pt(), event_weight);
                        h4mu_mu_eta_seen->Fill(vec_lepton.Eta(), event_weight);

                        //plotting acceptance against different variables for muons in only 4mu events seen by detector
                        e4mu_eta->FillWeighted(muon_cut, event_weight, vec_lepton.Eta());
                        e4mu_Et->FillWeighted(muon_cut, event_weight, TMath::Sqrt(vec_lepton.Pt() * vec_lepton.Pt() + vec_lepton.M() * vec_lepton.M()));
                        e4mu_pT->FillWeighted(muon_cut, event_weight, vec_lepton.Pt());

                        //plotting scatter graph for muons in 4mu events that are seen by the detector
                        h_4mu_pT_eta->Fill(vec_lepton.Eta(), vec_lepton.Pt());
                    }

                    if(abs(lep->PID) == 12){
                        h4mu_nu_pT_seen->Fill(vec_lepton.Pt(), event_weight);
                        h4mu_nu_eta_seen->Fill(vec_lepton.Eta(), event_weight); 
                    }
                }

                if(event4e_seen){
                    if(abs(lep->PID) == 11){
                        all_electrons_seen.push_back(lep);
                    }
                }

                if(event2e2mu_seen){
                    if(abs(lep->PID) == 11){
                        all_electrons_seen.push_back(lep);
                    }
                    if(abs(lep->PID) == 13){
                        all_muons_seen.push_back(lep);
                    }
                }
            }

            if(Debug){
                std::cout << " Missing Energy Four Vector: pT = " << missing_energy_vector.Pt() << " Et = " << missing_energy_vector.Et() 
                << " E = " << missing_energy_vector.E() << " eta = " << missing_energy_vector.Eta() << " phi = " << missing_energy_vector.Phi() 
                << " mass = " << missing_energy_vector.M() << " px = " << missing_energy_vector.Px() << " py = " << missing_energy_vector.Py() 
                << " pz = " << missing_energy_vector.Pz() << std::endl;
            }

            //------------------------------------------------------------------
            // Loop For 4l Events That Are Seen By The Detector
            //------------------------------------------------------------------
            if(event4mu_seen || event4e_seen || event2e2mu_seen){
                if(Debug){
                    std::cout << " muons seen = " << all_muons_seen.size() << std::endl;
                    std::cout << " electrons seen = " << all_electrons_seen.size() << std::endl;
                }

                missing_energy_vector = -missing_energy_vector;
                //missing energy plots
                h_ME_nu_pT -> Fill(vec_neutrino.Pt(), missing_energy_vector.Pt());
                h_ME_nu_eta -> Fill(vec_neutrino.Eta(), missing_energy_vector.Eta());
                h_ME_nu_phi -> Fill(vec_neutrino.Phi(), missing_energy_vector.Phi());

                //---------------------------------------------------------------
                // Kinematic Reconstruction
                //---------------------------------------------------------------

                std::vector<double> electron_output = Electron_Reconstruction(vec_neutrino);
                double Q2_e = electron_output[0];
                double x_e = electron_output[1];
                double y_e = electron_output[2];

                if(Debug){
                    std::cout << "Electron Reconstruction Method: Q_squared = " << Q2_e << " x = " << x_e << " y = " << y_e << std::endl;
                }

                x_Qsquared_electron -> Fill(x_e, Q2_e);
                h_logQsquared_electron -> Fill(TMath::Log10(Q2_e), event_weight);
                h_logx_electron -> Fill(TMath::Log10(x_e), event_weight);
                h_logy_electron -> Fill(TMath::Log10(y_e), event_weight);
                
                std::vector<double> hadron_output = Hadron_Reconstruction(missing_energy_vector);
                double Q2_h = hadron_output[0];
                double x_h = hadron_output[1];
                double y_h = hadron_output[2];

                if(Debug){
                    std::cout << "Hadron Reconstruction Method: Q_squared = " << Q2_h << " x = " << x_h << " y = " << y_h << std::endl;
                }

                x_Qsquared_hadron -> Fill(x_h, Q2_h);
                log_Qsquared_plot -> Fill(TMath::Log10(Q2_e), TMath::Log10(Q2_h));
                log_x_plot -> Fill(TMath::Log10(x_e), TMath::Log10(x_h));
                log_y_plot -> Fill(TMath::Log10(y_e), TMath::Log10(y_h));
                h_logQsquared_hadron -> Fill(TMath::Log10(Q2_h), event_weight);
                h_logx_hadron -> Fill(TMath::Log10(x_h), event_weight);
                h_logy_hadron -> Fill(TMath::Log10(y_h), event_weight);

                //---------------------------------------------------------------
                // Higgs and ZZ* Mass Reconstruction
                //---------------------------------------------------------------
                std::vector<GenParticle*> particles1;
                std::vector<GenParticle*> antiparticles1;
                std::vector<TLorentzVector> particles_vector;
                std::vector<TLorentzVector> antiparticles_vector;

                Particle_Antiparticle_Sorter(all_muons_seen, all_electrons_seen, particles1, antiparticles1);
                
                particles_vector = Make_Lorentz_Vector(particles1);
                antiparticles_vector = Make_Lorentz_Vector(antiparticles1);

                if(Debug) std::cout << " particles 1 size = " << particles1.size() << std::endl;
                if(Debug) std::cout << " antiparticles 1 size = " << antiparticles1.size() << std::endl;

                std::vector<double> output = Mass_Reconstruction(all_muons_seen, all_electrons_seen, particles_vector, antiparticles_vector);

                double m_4l = output[0];
                double Z_onshell = output[1];
                double Z_offshell = output[2];

                Zstar_cut_flags = Check_Cuts(Zstar_cut_values, Z_offshell, Zstar_cut_flags);
                Z_cut_flags = Check_Cuts(Z_cut_values, Z_onshell, Z_cut_flags);

                h_Higgs_reco->Fill(m_4l, event_weight);
                h_ZZ_mass_reco->Fill(Z_onshell, event_weight);
                h_ZZ_mass_reco->Fill(Z_offshell, event_weight);
                h_Z_reco->Fill(Z_onshell, event_weight);
                h_Zstar_reco->Fill(Z_offshell, event_weight);

                //---------------------------------------------------------------
                // Histograms for Significance and S/B Plots
                //---------------------------------------------------------------
                //fills each histogram with the reconstructed Higgs mass if the flag for that event is true
                Fill_Histogram(h_pT_cuts, pT_cut_flags, m_4l);    
                Fill_Histogram(h_Zstar_cuts, Zstar_cut_flags, m_4l);
                Fill_Histogram(h_Z_cuts, Z_cut_flags, m_4l);

                std::tuple<std::vector<TLorentzVector>, std::vector<int>> electron_smear;
                std::tuple<std::vector<TLorentzVector>, std::vector<int>> muon_smear;

                electron_smear = Smear(all_electrons_seen);
                muon_smear = Smear(all_muons_seen);

                std::tuple<std::vector<TLorentzVector>, std::vector<TLorentzVector>> particles_antiparticles_smeared;
                std::vector<TLorentzVector> particles_smeared;
                std::vector<TLorentzVector> antiparticles_smeared;
                
                particles_antiparticles_smeared = Smeared_Particle_Antiparticle_Sorter(electron_smear, muon_smear);
                particles_smeared = std::get<0>(particles_antiparticles_smeared);
                antiparticles_smeared = std::get<1>(particles_antiparticles_smeared);
                
                if(Debug) std::cout << "particles/antiparticles smeared size: " << particles_smeared.size() << antiparticles_smeared.size() << std::endl;

                std::vector<double> output1 = Mass_Reconstruction(all_muons_seen, all_electrons_seen, particles_smeared, antiparticles_smeared);

                double m_4l_smeared = output1[0];
                double Z_onshell_smeared = output1[1];
                double Z_offshell_smeared = output1[2];

                if(Debug) std::cout << "HIGGS MASS SMEARED: " << m_4l_smeared << std::endl; 

                h_Higgs_reco_smeared->Fill(m_4l_smeared, event_weight);
                h_Z_reco_smeared->Fill(Z_onshell_smeared, event_weight);
                h_Zstar_reco_smeared->Fill(Z_offshell_smeared, event_weight);

                // if(Debug){
                //     std::cout << " electron smear: px = " << electron_smear[0].Px() << " py = " << electron_smear[0].Py() << 
                //     " pz = " << electron_smear[0].Pz() << " E = " << electron_smear[0].E() << std::endl;
                // }

                // if(Debug){
                //     std::cout << " muon smear: px = " << muon_smear[0].Px() << " py = " << muon_smear[0].Py() << 
                //     " pz = " << muon_smear[0].Pz() << " E = " << muon_smear[0].E() << std::endl;
                // }
            }         
        }
        
        //------------------------------------------------------------------
        // Z/W/H Boson Loop
        //------------------------------------------------------------------
        std::vector<TLorentzVector> list_Zboson;

        for(int i = 0; i <bTruthWZ->GetEntriesFast(); ++i){
            
            GenParticle* boson = (GenParticle*) bTruthWZ->At(i);

            if(Debug) std::cout << "Boson " << i << " pT = " << boson->PT << " eta = " << boson->Eta << " phi = " 
            << boson->Phi << " mass = " << boson->Mass << std::endl;

            TLorentzVector vec_boson;
            vec_boson.SetPtEtaPhiM(boson->PT,boson->Eta,boson->Phi,boson->Mass);

            //plotting histograms for Z bosons in all events
            if(abs(boson->PID) == 23){
		        h_Z_Pt->Fill(vec_boson.Pt(), event_weight);
                h_Z_eta->Fill(vec_boson.Eta(), event_weight);
                h_Z_Et->Fill(TMath::Sqrt(vec_boson.Pt() * vec_boson.Pt() + vec_boson.M() * vec_boson.M()), event_weight);
                list_Zboson.push_back(vec_boson);
            }

            //looking for Higgs bosons
            if(boson->Mass == 125){
		        h_Higgs_pT->Fill(vec_boson.Pt(), event_weight);
                h_Higgs_eta->Fill(vec_boson.Eta(), event_weight);
                h_Higgs_Et->Fill(TMath::Sqrt(vec_boson.Pt() * vec_boson.Pt() + vec_boson.M() * vec_boson.M()), event_weight);
                e_H_eta->FillWeighted(muon_cut, event_weight, vec_boson.Eta());
                e_H_Et->FillWeighted(muon_cut, event_weight, TMath::Sqrt(vec_boson.Pt() * vec_boson.Pt() + vec_boson.M() * vec_boson.M()));
                e_H_pT->FillWeighted(muon_cut, event_weight, vec_boson.Pt());
                }

            if(event4mu_seen){
                if(abs(boson->PID) == 23){
                    h4mu_Z_pT_seen->Fill(vec_boson.Pt(), event_weight);
                    h4mu_Z_eta_seen->Fill(vec_boson.Eta(), event_weight);
                }

                if(boson->Mass == 125){
                    h4mu_Higgs_pT_seen->Fill(vec_boson.Pt(), event_weight);
                    h4mu_Higgs_eta_seen->Fill(vec_boson.Eta(), event_weight);
                }
            }
        }
        
        //Z/W/H Boson Loop End

        if(list_Zboson.size() > 1 ){

            TLorentzVector Higgs = list_Zboson.at(0) + list_Zboson.at(1);
                h_ZZ_Mass->Fill(Higgs.M(), event_weight);
        }
    } 
}

     // Loop over all events end

//-----------------------------------------------------------------------
// Kinematic Reconstruction
//-----------------------------------------------------------------------
std::vector<double> Electron_Reconstruction(TLorentzVector nu){
    std::vector<double> output;

    double root_s = 1700; // root(s) = 1.7TeV -> given in GeV, found from CDR update
    double e_E = 60; //in GeV, incoming electron energy
    double nu_E = abs(nu.E()); //scattered neutrino energy
    double angle = TMath::Pi() - nu.Theta();

    double Q2 = 2 * e_E * nu_E * (1 + TMath::Cos(angle));
    double y = 1 - (nu_E/e_E) * TMath::Power(TMath::Sin((angle)/2), 2);
    double x = Q2/(TMath::Power(root_s, 2) * y);

    output = {Q2, x, y};
    return output;
}

std::vector<double> Hadron_Reconstruction(TLorentzVector missing_energy_vector){
    std::vector<double> output;

    double root_s = 1700; // root(s) = 1.7TeV -> given in GeV, found from CDR update
    double incoming_e_E = 60; //in GeV, incoming electron energy
    double E = abs(missing_energy_vector.E());
    double pz = abs(missing_energy_vector.Pz());
    double pt = abs(missing_energy_vector.Pt());

    double sigma = E - pz;
    double y = sigma / (2 * incoming_e_E); //inelasticity for hadron reconstruction method
    double Q2 = (TMath::Power(pt, 2)) / (1 - y); //Q^2 for hadron reconstruction method
    double x = Q2/(TMath::Power(root_s, 2) * y); //Bjorken x for hadron reconstruction method

    output = {Q2, x, y};
    return output;
}

//-----------------------------------------------------------------------
// Defines Histograms for Significance and S/B Plots
//-----------------------------------------------------------------------
std::vector<TH1D*> Define_Histograms(TString hist_name, int n_cuts){
	TString suffix = "_cut";
    std::vector<TH1D*> h_varycuts;
	for(int n = 0; n < n_cuts; ++n){	
		TH1D * h_new = new TH1D(hist_name + suffix + std::to_string(n), " ; m_{4l} [Gev]; Events", 50, 50, 250);
		h_varycuts.push_back(h_new);
	}
    return h_varycuts;
}

//-----------------------------------------------------------------------
// Initialises Flags To Indicate Whether Each Event Passes Each Cut
//-----------------------------------------------------------------------
std::vector<bool> Initialise_Flags(int n_cuts){
    std::vector<bool> cut_flags;
    for(int n = 0; n < n_cuts; ++n){
        cut_flags.push_back(true);
    }
    return cut_flags;
}

//-----------------------------------------------------------------------
// Checks Whether Each Event Passes Each Cut and Changes Flag Accordingly
//-----------------------------------------------------------------------
std::vector<bool> Check_Cuts(std::vector<double> cut_values, double lepton_property, std::vector<bool> cut_flags){
    for(int i = 0; i < cut_values.size(); ++i){
        if(lepton_property < cut_values[i]){
            cut_flags[i] = false;
        }
    }
    return cut_flags;
}

//-----------------------------------------------------------------------
// Fills A Histogram For Each Cut Value With Events That Pass Cut
//-----------------------------------------------------------------------
void Fill_Histogram(std::vector<TH1D*> h_varycuts, std::vector<bool> cut_flags, double reco_Higgs){
    for(int i = 0; i < cut_flags.size(); ++i){
        if(cut_flags[i] == true){
            h_varycuts[i] -> Fill(reco_Higgs);
        }
    }
}

//-----------------------------------------------------------------------
// Writes n_cuts Number of Histograms
//-----------------------------------------------------------------------
void Write_Histogram(std::vector<TH1D*> h_varycuts){
    for(int i = 0; i < h_varycuts.size(); ++i){
        h_varycuts[i]->Write();
    }
}

//-----------------------------------------------------------------------
// Higgs and ZZ* Mass Reconstruction
//-----------------------------------------------------------------------
std::vector<TLorentzVector> Make_Lorentz_Vector(std::vector<GenParticle*> particles){
    TLorentzVector temp_vector;
    std::vector<TLorentzVector> output;

    for(int i = 0; i < particles.size(); ++i){
        temp_vector.SetPtEtaPhiM(particles[i]->PT, particles[i]->Eta, particles[i]->Phi, particles[i]->Mass);
        output.push_back(temp_vector);
    }

    return output;
}

void Particle_Antiparticle_Sorter(std::vector<GenParticle*> all_muons_seen, std::vector<GenParticle*> all_electrons_seen, std::vector<GenParticle*> &particles, std::vector<GenParticle*> &antiparticles){
    TLorentzVector temp_vector;

    if(all_muons_seen.size() == 4){
        for(int i = 0; i < all_muons_seen.size(); ++i){
            //temp_vector.SetPtEtaPhiM(all_muons_seen[i]->PT, all_muons_seen[i]->Eta, all_muons_seen[i]->Phi, all_muons_seen[i]->Mass);
            if(all_muons_seen[i]->PID == 13){
                particles.push_back(all_muons_seen[i]);
            }
            else if(all_muons_seen[i]->PID == -13){
                antiparticles.push_back(all_muons_seen[i]);
            }
        }  
    }

    if(all_electrons_seen.size() == 4){
        for(int i = 0; i < all_electrons_seen.size(); ++i){
            //temp_vector.SetPtEtaPhiM(all_electrons_seen[i]->PT, all_electrons_seen[i]->Eta, all_electrons_seen[i]->Phi, all_electrons_seen[i]->Mass);
            if(all_electrons_seen[i]->PID == 11){
                particles.push_back(all_electrons_seen[i]);
            }
            else if(all_electrons_seen[i]->PID == -11){
                antiparticles.push_back(all_electrons_seen[i]);
            }
        }  
    }

    if(all_electrons_seen.size() == 2 && all_muons_seen.size() == 2){
        for(int i = 0; i < all_electrons_seen.size(); ++i){
            //temp_vector.SetPtEtaPhiM(all_electrons_seen[i]->PT, all_electrons_seen[i]->Eta, all_electrons_seen[i]->Phi, all_electrons_seen[i]->Mass);
            if(all_electrons_seen[i]->PID == 11){
                particles.push_back(all_electrons_seen[i]);
            }
            else if(all_electrons_seen[i]->PID == -11){
                antiparticles.push_back(all_electrons_seen[i]);
            }
        }  
        for(int i = 0; i < all_muons_seen.size(); ++i){
            //temp_vector.SetPtEtaPhiM(all_muons_seen[i]->PT, all_muons_seen[i]->Eta, all_muons_seen[i]->Phi, all_muons_seen[i]->Mass);
            if(all_muons_seen[i]->PID == 13){
                particles.push_back(all_muons_seen[i]);
            }
            else if(all_muons_seen[i]->PID == -13){
                antiparticles.push_back(all_muons_seen[i]);
            }
        }
    }
}

std::tuple<std::vector<TLorentzVector>, std::vector<TLorentzVector>> Smeared_Particle_Antiparticle_Sorter(std::tuple<std::vector<TLorentzVector>, std::vector<int>> electron_smear, std::tuple<std::vector<TLorentzVector>, std::vector<int>> muon_smear){
    std::vector<TLorentzVector> electron_smear0;
    std::vector<int> electron_smear1;
    std::vector<TLorentzVector> muon_smear0;
    std::vector<int> muon_smear1;
    std::vector<TLorentzVector> particles_smeared;
    std::vector<TLorentzVector> antiparticles_smeared;

    electron_smear0 = std::get<0>(electron_smear);
    electron_smear1 = std::get<1>(electron_smear);
    muon_smear0 = std::get<0>(muon_smear);
    muon_smear1 = std::get<1>(muon_smear);

    for(int i = 0; i < electron_smear1.size(); ++i){
        if(electron_smear1[i] == 11){
            particles_smeared.push_back(electron_smear0[i]);
        }

        if(electron_smear1[i] == -11){
            antiparticles_smeared.push_back(electron_smear0[i]);
        }
    }

    for(int i = 0; i < muon_smear1.size(); ++i){
        if(muon_smear1[i] == 13){
            particles_smeared.push_back(muon_smear0[i]);
        }

        if(muon_smear1[i] == -13){
            antiparticles_smeared.push_back(muon_smear0[i]);
        }
    }

    return make_tuple(particles_smeared, antiparticles_smeared);
}

std::vector<double> Mass_Reconstruction(std::vector<GenParticle*> all_muons_seen, std::vector<GenParticle*> all_electrons_seen, std::vector<TLorentzVector> particles, std::vector<TLorentzVector> antiparticles){
    std::vector<double> output;
    TLorentzVector reco_Higgs;
    reco_Higgs.SetPtEtaPhiM(0,0,0,0);
    std::vector<TLorentzVector> recoZ; //makes an array which will contain all Lorentz vectors of reconstructed Z bosons from lepton pairs
    std::vector<double> recoZmass; //makes an array which will contain all potential reconstructed Z masses from lepton pairs
    std::vector<double> massdiff; //makes an array which will contain the difference between each reconstructed Z mass and the known Z mass
    double Zmass = 91.1876; //actual Z boson mass in GeV

    double Z_onshell;
    double Z_offshell;

    if(all_electrons_seen.size() == 2 && all_muons_seen.size() == 2){
        recoZ.push_back(particles[0]+antiparticles[0]);
        recoZ.push_back(particles[1]+antiparticles[1]);

        for(int j = 0; j < recoZ.size(); ++j){
            recoZmass.push_back(recoZ[j].M()); //we now have a list of all possible reconstructed Z masses
        }

        if(Debug){
            std::cout << " Z masses size " << recoZmass.size() << std::endl;
            std::cout << " Z masses: " << recoZmass[0] << "," << recoZmass[1] << std::endl;
        }

        for(int k = 0; k < recoZmass.size(); ++k){
            double diff = abs(recoZmass[k] - Zmass); //creates a list of the difference between the reconstructed Z mass and the known Z mass
            massdiff.push_back(diff);
        } 

        if(Debug){
            std::cout << " Mass Diff: " << massdiff[0] << "," << massdiff[1] << std::endl;
        }

        //the on-shell Z boson is the one that is closest to the known Z mass
        //from the lepton pair combinations above, once you know the on-shell Z, you can determine which is the off-shell Z
        double min = massdiff[0];
        for(int l = 0; l<massdiff.size(); ++l){ //finds the smallest mass difference
            if(massdiff[l] < min){
                min = massdiff[l];
            }
        }

        //these statements are all defining the on-shell and off-shell Z bosons based on the lepton pairings
        if(min == massdiff[0]){
            Z_onshell = recoZmass[0];
            Z_offshell = recoZmass[1];
        }

        if(min == massdiff[1]){
            Z_onshell = recoZmass[1];
            Z_offshell = recoZmass[0];
        }

        if(Debug){
        std::cout << " Smallest Mass Diff: " << min << std::endl;
        std::cout << " Z mass (on shell): " << Z_onshell << " Z mass (off shell): " << Z_offshell << std::endl;
    }
        
    } else if(all_electrons_seen.size() == 4 || all_muons_seen.size() == 4){
        recoZ.push_back(particles[0]+antiparticles[0]);
        recoZ.push_back(particles[0]+antiparticles[1]);
        recoZ.push_back(particles[1]+antiparticles[0]);
        recoZ.push_back(particles[1]+antiparticles[1]);

        for(int j = 0; j < recoZ.size(); ++j){
            recoZmass.push_back(recoZ[j].M()); //we now have a list of all possible reconstructed Z masses
        }

        if(Debug){
            std::cout << " Z masses size " << recoZmass.size() << std::endl;
            std::cout << " Z masses: " << recoZmass[0] << "," << recoZmass[1] << "," << recoZmass[2] << "," << recoZmass[3] << std::endl;
        }

        for(int k = 0; k < recoZmass.size(); ++k){
            double diff = abs(recoZmass[k] - Zmass); //creates a list of the difference between the reconstructed Z mass and the known Z mass
            massdiff.push_back(diff);
        } 

        if(Debug){
            std::cout << " Mass Diff: " << massdiff[0] << "," << massdiff[1] << "," << massdiff[2] << "," << massdiff[3] << std::endl;
        }

        //the on-shell Z boson is the one that is closest to the known Z mass
        //from the lepton pair combinations above, once you know the on-shell Z, you can determine which is the off-shell Z
        double min = massdiff[0];
        for(int l = 0; l<massdiff.size(); ++l){ //finds the smallest mass difference
            if(massdiff[l] < min){
                min = massdiff[l];
            }
        }

        //these statements are all defining the on-shell and off-shell Z bosons based on the lepton pairings
        if(min == massdiff[0]){
            Z_onshell = recoZmass[0];
            Z_offshell = recoZmass[3];
        }

        if(min == massdiff[1]){
            Z_onshell = recoZmass[1];
            Z_offshell = recoZmass[2];
        }

        if(min == massdiff[2]){
            Z_onshell = recoZmass[2];
            Z_offshell = recoZmass[1];
        }

        if(min == massdiff[3]){
            Z_onshell = recoZmass[3];
            Z_offshell = recoZmass[0];
        }

        if(Debug){
            std::cout << " Smallest Mass Diff: " << min << std::endl;
            std::cout << " Z mass (on shell): " << Z_onshell << " Z mass (off shell): " << Z_offshell << std::endl;
        }
    }

    reco_Higgs = particles[0] + particles[1] + antiparticles[0] + antiparticles[1];

    if(Debug){
        std::cout << " Higgs Mass: " << reco_Higgs.M() << std::endl;
        std::cout << " particles.size = " << particles.size() << std::endl;
        std::cout << " antiparticles.size = " << antiparticles.size() << std::endl;
    }

    output = {reco_Higgs.M(), Z_onshell, Z_offshell};
    return output;
}

//------------------------------------------------------------------
// Smearing the Energy and Momentum for Electrons
//------------------------------------------------------------------
std::tuple<std::vector<TLorentzVector>, std::vector<int>> Smear(std::vector<GenParticle*> particles){
    std::tuple<std::vector<TLorentzVector>, std::vector<int>> output;
    std::vector<TLorentzVector> vectors;
    std::vector<int> PIDs;
    std::vector<double> random_num_list;

    TRandom3 * generate_random_num = new TRandom3();
    generate_random_num -> SetSeed(0);

    double a = 12.4/100; //found from CDR update
    double b = 1.9/100;

    for(int i = 0; i < particles.size(); ++i){
        double E = particles[i]->E;
        double px = particles[i]->Px;
        double py = particles[i]->Py;
        double pz = particles[i]->Pz;
        double random_num = generate_random_num -> Gaus(0, 1);

        if(abs(particles[i]->PID) == 11){
            //random_num_list.push_back(random_num);
            double E_resolution = E * TMath::Sqrt(TMath::Power((a/TMath::Sqrt(E)), 2) + TMath::Power(b, 2));
            //double E_resolution = 0;
            double E_smear = E + E_resolution * random_num;
            double px_smear = px + E_resolution * random_num;
            double py_smear = py + E_resolution * random_num;
            double pz_smear = pz + E_resolution * random_num;

            TLorentzVector electron_vector;
            electron_vector.SetPxPyPzE(px_smear, py_smear, pz_smear, E_smear);
            vectors.push_back(electron_vector);
            PIDs.push_back(particles[i]->PID);
        }

        if(abs(particles[i]->PID) == 13){
            double p_resolution = 0.02; //the momentum resolution of the inner tracker is 1-2%, I have used 2% here as this is a "worse case scenario" of the detector effects
            //double p_resolution = 0; 
            double E_smear = E + p_resolution * E * random_num;
            double px_smear = px + p_resolution * px * random_num;
            double py_smear = py + p_resolution * py * random_num;
            double pz_smear = pz + p_resolution * pz * random_num;

            TLorentzVector muon_vector;
            muon_vector.SetPxPyPzE(px_smear, py_smear, pz_smear, E_smear);
            vectors.push_back(muon_vector);
            PIDs.push_back(particles[i]->PID);
        }
    }
    
    return make_tuple(vectors, PIDs);
}

//------------------------------------------------------------------
// Smearing the Energy and Momentum for Muons
//------------------------------------------------------------------

// void PaintBin (TH1D * histogram, Int_t bin, Int_t color){
//    printf("%d %d %d\n", bin, color, histogram->GetBinContent(bin));
//    TBox * box = new TBox(histogram->GetBinLowEdge(bin),
//                          histogram->GetMinimum(),
//                          histogram->GetBinWidth(bin)+histogram->GetBinLowEdge(bin),
//                          histogram->GetBinContent(bin));
//    box->SetFillColor(color);
//    box->Draw();
// }