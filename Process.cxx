
#include "Process.h"
#include "utilities.cxx"

bool Debug = false;
//const TString signal_inputfile = "/disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_H4l.Delphes.root";
//const TString bkgd_inputfile = "/disk/moose/general/asc/Y4_LHeC_2021/Samples_Py8303_Delphes3.4.3pre06/LHeC_CC_Zll.100k.v2.Delphes.root";

int main(int argc, char* argv[]) {

    // Input Delphes File
 
    const TString InputFile = argv[1];
    const TString OutputFileName = argv[2];

    //const bool signal = argv[1] == "signal.root";

    //const TString InputFile = signal ? signal_inputfile : bkgd_inputfile;

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

    //this will be a histogram where the first bin is the total number of events, the second bin is the number of 4mu
    //events and the third bin is the number of 4mu events that pass the cuts
    h_EventCount = new TH1D("h_EventCount",";""; Number of Events",7,0,7);
    TAxis * xAxis = h_EventCount->GetXaxis();
    xAxis->SetBinLabel(1, "Total N Events");
    xAxis->SetBinLabel(2, "N 4mu Events");
    xAxis->SetBinLabel(3, "N 4mu Events Seen by Detector");
    xAxis->SetBinLabel(4, "N 4e Events");
    xAxis->SetBinLabel(5, "N 2e2mu Events");
    xAxis->SetBinLabel(6, "N 4l Events");
    xAxis->SetBinLabel(7, "N 4l Events Seen by Detector");


    h_WeightCount = new TH1D("h_WeightCount","",1,0,1);

    h_ZZ_Mass = new TH1D("h_ZZ_Mass","; ZZ Mass [GeV]; Events / 2 GeV",125,0.0,250.0);

    //transverse momentum histograms
    h_Lepton_Pt = new TH1D("h_Lepton_Pt","; Charged Lepton p_{T} [GeV]; Events / 5 GeV",200,0.0,300.0);
    h_Z_Pt = new TH1D("h_Z_Pt","; Z p_{T} [GeV]; Events / 5 GeV",200,0.0,300.0);
    h_Jet_Pt = new TH1D("h_Jet_Pt","; Jet p_{T} [GeV]; Events / 5 GeV",200,0.0,300.0);
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
        //comparing electron and hadron reconstruction methods
    log_Qsquared_plot = new TH2D("log_Qsquared_plot", "log_{10}(Q^{2}) (Electron) against log_{10}(Q^{2}) (Hadron); log_{10}(Q^{2}_{h}); log_{10}(Q^{2}_{e})", 100, 0, 7, 100, 0, 7);
    log_x_plot = new TH2D("log_x_plot", "log_{10}(x) (Electron) against log_{10}(x) (Hadron); log(x_{h}); log_{10}(x_{e})", 100, -5, 0, 100, -5, 0);
    log_y_plot = new TH2D("log_y_plot", "log_{10}(y) (Electron) against log_{10}(y) (Hadron); log(y_{h}); log_{10}(y_{e})", 100, -2, 1, 100, -2, 1);

    //mass reconstruction
    h_Higgs_reco = new TH1D("h_Higgs_reco", "; Reconstructed Higgs Boson Mass [GeV]; Events ", 150, 0.0, 130.0);
    h_ZZ_mass_reco = new TH1D("h_ZZ_mass_reco", "; Reconstructed ZZ* Mass [GeV]; Events ", 150, 0.0, 130.0);
    h_Z_reco = new TH1D("h_Z_reco", "; Reconstructed Z Boson Mass [GeV]; Events ", 150, 0.0, 130.0);
    h_Zstar_reco = new TH1D("h_Zstar_reco", "; Reconstructed Z* Boson Mass [GeV]; Events ", 150, 0.0, 130.0);

    //event weight plot
    h_eventweight = new TH1D("h_eventweight", "; Event Weight ; Events ", 100, 0.000100, 0.000120);

    //mass reconstruction plots with cuts implemented
    // h_Zstar_cuts = new TH1D("h_Zstar_cuts", "; Cut on Maximum Z* Mass [GeV] ; Events ", 150, 0.0, 100);
    

    // Run the selection
    Process(reader); //removed signal from the arguments

    std::cout << "Total Number of Events: " << h_EventCount->GetBinContent(1) << std::endl;
    std::cout << "Number of 4mu Events: " << h_EventCount->GetBinContent(2) << std::endl;
    std::cout << "Number of 4mu Events Seen by Detector: " << h_EventCount->GetBinContent(3) << std::endl; 
    std::cout << "Number of 4e Events: " << h_EventCount->GetBinContent(4) << std::endl;
    std::cout << "Number of 2e2mu Events: " << h_EventCount->GetBinContent(5) << std::endl;
    std::cout << "Number of 4l Events: " << h_EventCount->GetBinContent(6) << std::endl;
    std::cout << "Number of 4l Events Seen by Detector: " << h_EventCount->GetBinContent(7) << std::endl;

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

    OutputFile->cd("Transverse Energy"); //maybe recreate these with total energy?

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
    h_Higgs_reco->GetXaxis()->SetRangeUser(0, 130);
    h_Higgs_reco->SetLineColor(kRed);
    h_Higgs_reco->SetStats(kFALSE);
    h_Higgs_reco->Draw("hist E2");

    h_Z_reco->SetLineColor(kGreen);
    h_Z_reco->Draw("hist same E2");

    h_Zstar_reco->SetLineColor(kBlue);
    h_Zstar_reco->Draw("hist same E2");

    legend = new TLegend(0.8,0.8,0.9,0.9);
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
    // for(int i = 0; i < h_pT_cuts.size(); ++i){
    //     h_pT_cuts[i]->Write();
    // }

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


void Process(ExRootTreeReader * treeReader) { //removed bool signal from the arguments

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
    
    //---------------------------------------------------------
    //calculating scale factors for signal and background files
    //---------------------------------------------------------
    double sigma_signal = 1.34e-17; //in barns
    double sigma_bkgd = 8.9e-18; //in barns
    double no_events = 100000;
    double luminosity_signal = no_events/sigma_signal;
    double luminosity_bkgd = no_events/sigma_bkgd;
    double luminosity_LHeC = 1000e15; //1000 inverse femtobarns
    double sf_signal = luminosity_LHeC/luminosity_signal; //scale factor
    double sf_bkgd = luminosity_LHeC/luminosity_bkgd;

    std::cout << "sf signal = " << sf_signal << " sf bkgd = " << sf_bkgd << std::endl;

    std::vector<bool> pT_cuts; //defines an array which will contain a flag to indicate whether each event passes each minimum pT cut

    h_pT_cuts = Define_Histograms("h_pT", 40);
    h_Zstar_cuts = Define_Histograms("h_Zstar", 40);

    std::vector<double> pT_cut_values = Define_Cut_Values(40, 10.0, 4.0);
    std::vector<double> Zstar_cut_values = Define_Cut_Values(40, 10.0, 10.0);
    
    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        //Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0);

        //if the argument = signal, it will use the signal scale factor. if not, it will use the background scale factor
        //const float Event_Weight = signal ? sf_signal : sf_bkgd;
        const float Event_Weight = sf_signal; //manually change this

        h_EventCount->Fill(0.5);
        h_WeightCount->Fill(0.5,Event_Weight);
        h_eventweight->Fill(Event_Weight);

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        bool event4mu = false; //flag to check if the event satisfies the 4mu subchannel
        bool event4mu_seen = true; //flag to check if the 4mu event can be seen by the detector
        bool event4e = false;
        bool event4e_seen = true;
        bool event2e2mu = false;
        bool event2e2mu_seen = true;

        //This is to calculate the missing energy, set it equal to zero at the start of each event and the vector 
        //sum will accumulate throughout.
        TLorentzVector Missing_Energy_Vector;
        Missing_Energy_Vector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);

        //------------------------------------------------------------------
        // Jet Loop
        //------------------------------------------------------------------

        bool JetCut = false; //flag to check if the jet passes the jet cut
        double eta_min_j = -5.0; //eta range of the hadronic calorimeter
        double eta_max_j = 5.5;
        double pT_min = 15; //in GeV

        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet * jet = (Jet*) bJet->At(i);

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << 
            " mass = " << jet->Mass << " flavour = " << jet->Flavor << std::endl;

            TLorentzVector Vec_Jet;
            Vec_Jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
            
            //defining the selection criteria for the jets
            if(
                eta_min_j<jet->Eta && 
                jet->Eta<eta_max_j && 
                jet->PT>pT_min){
                JetCut = true;
                }

            if(JetCut){
                h_Jet_Pt->Fill(Vec_Jet.Pt(), Event_Weight);
                h_Jet_eta->Fill(Vec_Jet.Eta(), Event_Weight);
                h_Jet_Et->Fill(TMath::Sqrt(Vec_Jet.Pt() * Vec_Jet.Pt() + Vec_Jet.M() * Vec_Jet.M()), Event_Weight);
            }

            TLorentzVector Vec_Lep;

            bool truejet = true; //flag to check if the jet is a true jet or not

            //lepton loop inside jet loop to calculate DeltaR for each jet/lepton combination
            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) { 

                GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

                Vec_Lep.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);

                double deltaR = TMath::Power(Vec_Jet.Phi() - Vec_Lep.Phi(), 2) + TMath::Power(Vec_Jet.Eta() - Vec_Lep.Eta(), 2);
                double jet_radius = 0.4; //defined in the DELPHES configuration

                //a true jet is defined as one where the DeltaR between the jet and the lepton is greater than the jet radius
                if(deltaR < jet_radius) truejet = false;

                // if(Debug) std::cout << "  DeltaR (Lepton " << i << ") = " << deltaR << std::endl;

            }

            if(Debug) std::cout << "  Is it a true jet? " << (truejet ? "yes" : "no") << std::endl;
        
            //jet loop again but for only true jets that pass the jet cuts
            if(truejet && JetCut){

                //adding true jet four vectors to the missing energy four vector
                Missing_Energy_Vector = Missing_Energy_Vector + Vec_Jet; 
                if(Debug) std::cout << "  Missing Energy Vector (Jets) " << Missing_Energy_Vector.Pt() << std::endl;

                h_trueJet_Pt->Fill(Vec_Jet.Pt(), Event_Weight);
                h_trueJet_eta->Fill(Vec_Jet.Eta(), Event_Weight);
                h_trueJet_Et->Fill(TMath::Sqrt(Vec_Jet.Pt() * Vec_Jet.Pt() + Vec_Jet.M() * Vec_Jet.M()), Event_Weight);
                
            }
            
        }
         // Jet Loop

        bool MuonCut = true; //flag to check if muons pass cuts
        double eta_min_mu = -4.6; //eta range of the inner tracker
        double eta_max_mu = 5.3;
        //momentum cut has been previously defined
        bool ElectronCut = true;
        double eta_min_e = -5.0;
        double eta_max_e = 5.5;

        std::vector<bool> pT_cut_flags = Initialise_Flags(40);
        std::vector<bool> Zstar_cut_flags = Initialise_Flags(40);
        
        // pT_cuts.clear(); //empties the vector once each event
        // for(int n = 0; n < h_pT_cuts.size(); ++n){
        //     pT_cuts.push_back(true); //set all the flags equal to true
        // }

        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------
        std::vector<GenParticle*> all_muons; //makes an array which will contain all muons
        std::vector<GenParticle*> all_electrons;
        
        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

            if(Debug){
                std::cout << "Lepton " << i << " PID = " << lep->PID << " pT = " << lep->PT << " eta = " << lep->Eta
                << " phi = " << lep->Phi << " mass = " << lep->Mass << " energy = " << lep->E << " momentum = " << 
                lep->P << std::endl;
            }

            if(abs(lep->PID) == 11){ 
                all_electrons.push_back(lep);
            }

            if(abs(lep->PID) == 13){ 
                all_muons.push_back(lep);
            }

            if(all_electrons.size() == 4){
                event4e = true;
            }

            if(all_muons.size() == 4){
                event4mu = true;
            }

            if(all_electrons.size() == 2 && all_muons.size() == 2){
                event2e2mu = true;
            }

            //defining the selection criteria for the muons
            if(abs(lep->PID) == 13){
                if(
                    lep->Eta > eta_min_mu &&
                    lep->Eta < eta_max_mu && 
                    lep->PT > pT_min) {
                        MuonCut = true;
                    } else{
                        MuonCut = false;
                    }
                
                if(Debug) std::cout << "  Does it pass the muon cut? " << (MuonCut ? "yes" : "no") << std::endl;

                if(!MuonCut){
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
                        ElectronCut = true;
                    } else{
                        ElectronCut = false;
                    }
                
                if(Debug) std::cout << "  Does it pass the electron cut? " << (ElectronCut ? "yes" : "no") << std::endl;

                if(!ElectronCut){
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
        
            TLorentzVector Vec_Lepton1;
            Vec_Lepton1.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
            
            // Look for electrons or muons
            if(abs(lep->PID) == 11 || abs(lep->PID) == 13 ){
                    h_Lepton_Pt->Fill( Vec_Lepton1.Pt(), Event_Weight );  
            }

            //plotting histograms for muons that pass the cuts
            if(abs(lep->PID) == 13 && MuonCut){
                h_mu_eta -> Fill(Vec_Lepton1.Eta(), Event_Weight);
                h_mu_Et -> Fill (TMath::Sqrt(Vec_Lepton1.Pt() * Vec_Lepton1.Pt() + Vec_Lepton1.M() * Vec_Lepton1.M()), Event_Weight);
                h_mu_pT -> Fill(Vec_Lepton1.Pt(), Event_Weight);
            }

            //plotting acceptance against different variables for muons in all events
            if(abs(lep->PID) == 13){
                e_eta->FillWeighted(MuonCut, Event_Weight, Vec_Lepton1.Eta());
                e_Et->FillWeighted(MuonCut, Event_Weight, TMath::Sqrt(Vec_Lepton1.Pt() * Vec_Lepton1.Pt() + Vec_Lepton1.M() * Vec_Lepton1.M()));
                e_pT->FillWeighted(MuonCut, Event_Weight, Vec_Lepton1.Pt());
                e_mu_pT_eta->FillWeighted(Vec_Lepton1.Pt(), Event_Weight, Vec_Lepton1.Eta());
            }

            //looking for neutrinos
            if(abs(lep->PID) == 12){
                h_nu_Et -> Fill(TMath::Sqrt(Vec_Lepton1.Pt() * Vec_Lepton1.Pt() + Vec_Lepton1.M() * Vec_Lepton1.M()), Event_Weight);
                h_nu_eta -> Fill(Vec_Lepton1.Eta(), Event_Weight);
            }

            //if any electrons are seen in the decay then it is not a 4mu decay mode and so will be ignored
            // if(abs(lep->PID) == 11){ 
            //     event4mu = false;
            //     all_electrons.push_back(lep);
            // }

            // if(abs(lep->PID) == 13){ 
            //     event4e = false;
            //     all_muons.push_back(lep);
            // }

        } // Lepton Loop End

        // if(all_electrons.size() != 2 && all_muons.size() != 2){
        //     event2e2mu = false;
        // }

        if(event4mu){
            h_EventCount -> Fill(1.5);
            h_EventCount -> Fill(5.5);
        }

        if(event4e){
            h_EventCount->Fill(3.5);
            h_EventCount->Fill(5.5);
        }

        if(event2e2mu){
            h_EventCount->Fill(4.5);
            h_EventCount->Fill(5.5);
        }

        if(event4mu_seen){
            h_EventCount -> Fill(2.5);
            h_EventCount -> Fill(6.5);
        }

        if(event4e_seen){
            h_EventCount -> Fill(6.5);
        }

        if(event2e2mu_seen){
            h_EventCount -> Fill(6.5);
        }
        
        
        //another lepton loop but only including the 4mu decays
        if(event4mu || event4e || event2e2mu){

            if(Debug){
                std::cout << "  Is it a 4mu event? " << (event4mu ? "yes" : "no") << std::endl;
                std::cout << "  Is it a 4e event? " << (event4e ? "yes" : "no") << std::endl;
                std::cout << "  Is it a 2e2mu event? " << (event2e2mu ? "yes" : "no") << std::endl;
                std::cout << "  Is it a 4mu event seen by the detector? " << (event4mu_seen ? "yes" : "no") << std::endl;
                std::cout << "  Is it a 4e event seen by the detector? " << (event4e_seen ? "yes" : "no") << std::endl;
                std::cout << "  Is it a 2e2mu event seen by the detector? " << (event2e2mu_seen ? "yes" : "no") << std::endl;
            }
            
            TLorentzVector Vec_Neutrino;
            std::vector<GenParticle*> all_muons_seen; //makes an array which will contain all muons that are seen by the detector

            for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

                GenParticle * lep = (GenParticle*) bTruthLepton->At(i);

                TLorentzVector Vec_Lepton;
                Vec_Lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);

                //adding all of the muon four vectors to the missing energy four vector
                if(lep->Status == 1 && abs(lep->PID) == 13){
                    Missing_Energy_Vector = Missing_Energy_Vector + Vec_Lepton;
                }

                if(lep->Status == 1 && abs(lep->PID) == 11){
                    Missing_Energy_Vector = Missing_Energy_Vector + Vec_Lepton;
                }

                if(event4mu){
                    //plotting histograms for muons in 4mu events
                    if(abs(lep->PID) == 13){
                        h4mu_mu_eta -> Fill(Vec_Lepton.Eta(), Event_Weight);
                        h4mu_mu_Et -> Fill (TMath::Sqrt(Vec_Lepton.Pt() * Vec_Lepton.Pt() + Vec_Lepton.M() * Vec_Lepton.M()), Event_Weight);
                        h4mu_mu_pT-> Fill(Vec_Lepton.Pt(), Event_Weight);
                        h_mu_pT_eta -> Fill(Vec_Lepton.Eta(), Vec_Lepton.Pt());
                    }

                    //plotting histograms for neutrinos in 4mu events
                    if(abs(lep->PID) == 12){
                        Vec_Neutrino = Vec_Lepton;
                        h4mu_nu_pT -> Fill(Vec_Lepton.Pt(), Event_Weight);
                        h4mu_nu_eta -> Fill(Vec_Lepton.Eta(), Event_Weight); 
                    }
                }
                
                //plots for 4mu events that are seen by the detector
                if(event4mu_seen){
                    if(abs(lep->PID) == 13){
                        all_muons_seen.push_back(lep);
                        h4mu_mu_pT_seen->Fill(Vec_Lepton.Pt(), Event_Weight);
                        h4mu_mu_eta_seen->Fill(Vec_Lepton.Eta(), Event_Weight);

                        //plotting acceptance against different variables for muons in only 4mu events seen by detector
                        e4mu_eta->FillWeighted(MuonCut, Event_Weight, Vec_Lepton.Eta());
                        e4mu_Et->FillWeighted(MuonCut, Event_Weight, TMath::Sqrt(Vec_Lepton.Pt() * Vec_Lepton.Pt() + Vec_Lepton.M() * Vec_Lepton.M()));
                        e4mu_pT->FillWeighted(MuonCut, Event_Weight, Vec_Lepton.Pt());

                        //plotting scatter graph for muons in 4mu events that are seen by the detector
                        h_4mu_pT_eta->Fill(Vec_Lepton.Eta(), Vec_Lepton.Pt());
                    }

                    if(abs(lep->PID) == 12){
                        h4mu_nu_pT_seen->Fill(Vec_Lepton.Pt(), Event_Weight);
                        h4mu_nu_eta_seen->Fill(Vec_Lepton.Eta(), Event_Weight); 
                    }

                }
            }

            if(Debug){
                std::cout << " Missing Energy Four Vector: pT = " << Missing_Energy_Vector.Pt() << " Et = " << Missing_Energy_Vector.Et() 
                << " E = " << Missing_Energy_Vector.E() << " eta = " << Missing_Energy_Vector.Eta() << " phi = " << Missing_Energy_Vector.Phi() 
                << " mass = " << Missing_Energy_Vector.M() << " px = " << Missing_Energy_Vector.Px() << " py = " << Missing_Energy_Vector.Py() 
                << " pz = " << Missing_Energy_Vector.Pz() << std::endl;
            }

            if(event4mu_seen){
                Missing_Energy_Vector = -Missing_Energy_Vector;
                //missing energy plots
                h_ME_nu_pT -> Fill(Vec_Neutrino.Pt(), Missing_Energy_Vector.Pt());
                h_ME_nu_eta -> Fill(Vec_Neutrino.Eta(), Missing_Energy_Vector.Eta());
                h_ME_nu_phi -> Fill(Vec_Neutrino.Phi(), Missing_Energy_Vector.Phi());

                std::vector<double> electron_output = Electron_Reconstruction(Vec_Neutrino);
                double Q2_e = electron_output[0];
                double x_e = electron_output[1];
                double y_e = electron_output[2];

                if(Debug){
                    std::cout << "Electron Reconstruction Method: Q_squared = " << Q2_e << " x = " << x_e << " y = " << y_e << std::endl;
                }

                x_Qsquared_electron -> Fill(x_e, Q2_e);
                h_logQsquared_electron -> Fill(TMath::Log10(Q2_e), Event_Weight);
                h_logx_electron -> Fill(TMath::Log10(x_e), Event_Weight);
                h_logy_electron -> Fill(TMath::Log10(y_e), Event_Weight);
                
                std::vector<double> hadron_output = Hadron_Reconstruction(Missing_Energy_Vector);
                double Q2_h = hadron_output[0];
                double x_h = hadron_output[1];
                double y_h = hadron_output[2];

                x_Qsquared_hadron -> Fill(x_h, Q2_h);
                log_Qsquared_plot -> Fill(TMath::Log10(Q2_e), TMath::Log10(Q2_h));
                log_x_plot -> Fill(TMath::Log10(x_e), TMath::Log10(x_h));
                log_y_plot -> Fill(TMath::Log10(y_e), TMath::Log10(y_h));
                h_logQsquared_hadron -> Fill(TMath::Log10(Q2_h), Event_Weight);
                h_logx_hadron -> Fill(TMath::Log10(x_h), Event_Weight);
                h_logy_hadron -> Fill(TMath::Log10(y_h), Event_Weight);

                //reconstructing Z mass from lepton pairs
                std::vector<TLorentzVector> muons; //makes an array which will contain the Lorentz vectors for muons
                TLorentzVector reco_Higgs;
                reco_Higgs.SetPtEtaPhiM(0,0,0,0);
                std::vector<TLorentzVector> antimuons; //makes an array which will contain the Lorentz vectors for antimuons
                std::vector<TLorentzVector> recoZ; //makes an array which will contain all Lorentz vectors of reconstructed Z bosons from lepton pairs
                TLorentzVector Temp_Vector;
                std::vector<double> recoZmass; //makes an array which will contain all potential reconstructed Z masses from lepton pairs
                std::vector<double> massdiff; //makes an array which will contain the difference between each reconstructed Z mass and the known Z mass
                    
                for(int i = 0; i < all_muons_seen.size(); ++i){

                    Temp_Vector.SetPtEtaPhiM(all_muons_seen[i]->PT, all_muons_seen[i]->Eta, all_muons_seen[i]->Phi, all_muons_seen[i]->Mass);
                    if(all_muons_seen[i]->PID==13){
                        muons.push_back(Temp_Vector);
                    }
                    else if(all_muons_seen[i]->PID==-13){
                        antimuons.push_back(Temp_Vector);
                    }
                }

                if(Debug){
                    std::cout << " Here 1 " << muons[0].M() << muons[1].M() << std::endl;
                }

                reco_Higgs = muons[0] + muons[1] + antimuons[0] + antimuons[1];

                if(Debug){
                    std::cout << " Here 2 " << std::endl;
                }

                recoZ.push_back(muons[0]+antimuons[0]);
                recoZ.push_back(muons[0]+antimuons[1]);
                recoZ.push_back(muons[1]+antimuons[0]);
                recoZ.push_back(muons[1]+antimuons[1]); //we now have a list of the Lorentz Vectors from every combination of lepton pairs

                if(Debug){
                    std::cout << " Here 3 " << std::endl;
                }

                for(int j = 0; j < recoZ.size(); ++j){
                    recoZmass.push_back(recoZ[j].M()); //we now have a list of all possible reconstructed Z masses
                } 

                if(Debug){
                    std::cout << " Z masses: " << recoZmass[0] << "," << recoZmass[1] << "," << recoZmass[2] << "," << recoZmass[3] << std::endl;
                }

                double Zmass = 91.1876; //actual Z boson mass in GeV

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

                double Z_onshell;
                double Z_offshell;

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
                    std::cout << " Higgs Mass: " << reco_Higgs.M() << std::endl;
                } 

                Zstar_cut_flags = Check_Cuts(Zstar_cut_values, Z_offshell, Zstar_cut_flags);

                h_Higgs_reco->Fill(reco_Higgs.M(), Event_Weight);
                h_ZZ_mass_reco->Fill(Z_onshell, Event_Weight);
                h_ZZ_mass_reco->Fill(Z_offshell, Event_Weight);
                h_Z_reco->Fill(Z_onshell, Event_Weight);
                h_Zstar_reco->Fill(Z_offshell, Event_Weight);

                //plots for significance
                    //fills each histogram with the reconstructed Higgs mass if the flag for that event is true
                Fill_Histogram(h_pT_cuts, pT_cut_flags, reco_Higgs.M());    
                Fill_Histogram(h_Zstar_cuts, Zstar_cut_flags, reco_Higgs.M());


                // for(int i = 0; i < pT_cuts.size(); ++i){
                //     if(pT_cuts[i] == true){
                //         h_pT_cuts[i] -> Fill(reco_Higgs.M());
                //     }
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

            TLorentzVector Vec_Boson;
            Vec_Boson.SetPtEtaPhiM(boson->PT,boson->Eta,boson->Phi,boson->Mass);

            //plotting histograms for Z bosons in all events
            if(abs(boson->PID) == 23){
		        h_Z_Pt->Fill(Vec_Boson.Pt(), Event_Weight);
                h_Z_eta->Fill(Vec_Boson.Eta(), Event_Weight);
                h_Z_Et->Fill(TMath::Sqrt(Vec_Boson.Pt() * Vec_Boson.Pt() + Vec_Boson.M() * Vec_Boson.M()), Event_Weight);
                list_Zboson.push_back(Vec_Boson);
            }

            //looking for Higgs bosons
            if(boson->Mass == 125){
		        h_Higgs_pT->Fill(Vec_Boson.Pt(), Event_Weight);
                h_Higgs_eta->Fill(Vec_Boson.Eta(), Event_Weight);
                h_Higgs_Et->Fill(TMath::Sqrt(Vec_Boson.Pt() * Vec_Boson.Pt() + Vec_Boson.M() * Vec_Boson.M()), Event_Weight);
                e_H_eta->FillWeighted(MuonCut, Event_Weight, Vec_Boson.Eta());
                e_H_Et->FillWeighted(MuonCut, Event_Weight, TMath::Sqrt(Vec_Boson.Pt() * Vec_Boson.Pt() + Vec_Boson.M() * Vec_Boson.M()));
                e_H_pT->FillWeighted(MuonCut, Event_Weight, Vec_Boson.Pt());
                }

            if(event4mu_seen){
                if(abs(boson->PID) == 23){
                    h4mu_Z_pT_seen->Fill(Vec_Boson.Pt(), Event_Weight);
                    h4mu_Z_eta_seen->Fill(Vec_Boson.Eta(), Event_Weight);
                }

                if(boson->Mass == 125){
                    h4mu_Higgs_pT_seen->Fill(Vec_Boson.Pt(), Event_Weight);
                    h4mu_Higgs_eta_seen->Fill(Vec_Boson.Eta(), Event_Weight);
                }
            }
        }
        
        //Z/W/H Boson Loop

        if(list_Zboson.size() > 1 ){

            TLorentzVector Higgs = list_Zboson.at(0) + list_Zboson.at(1);
                h_ZZ_Mass->Fill(Higgs.M(), Event_Weight);
        }
    }
    
}

     // Loop over all events

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

std::vector<double> Hadron_Reconstruction(TLorentzVector Missing_Energy_Vector){
    std::vector<double> output;

    double root_s = 1700; // root(s) = 1.7TeV -> given in GeV, found from CDR update
    double incoming_e_E = 60; //in GeV, incoming electron energy
    double E = abs(Missing_Energy_Vector.E());
    double pz = abs(Missing_Energy_Vector.Pz());
    double pt = abs(Missing_Energy_Vector.Pt());

    double Sigma = E - pz;
    double y = Sigma / (2 * incoming_e_E); //inelasticity for hadron reconstruction method
    double Q2 = (TMath::Power(pt, 2)) / (1 - y); //Q^2 for hadron reconstruction method
    double x = Q2/(TMath::Power(root_s, 2) * y); //Bjorken x for hadron reconstruction method

    output = {Q2, x, y};
    return output;
}

std::vector<TH1D*> Define_Histograms(TString hist_name, int n_cuts){
	TString suffix = "_cut";
    std::vector<TH1D*> h_varycuts;
	for(int n = 0; n < n_cuts; ++n){	
		TH1D * h_new = new TH1D(hist_name + suffix + std::to_string(n), " ; m_{4l} [Gev]; Events", 50, 50, 250);
		h_varycuts.push_back(h_new);
	}
    return h_varycuts;
}

std::vector<bool> Initialise_Flags(int n_cuts){
    std::vector<bool> cut_flags;
    for(int n = 0; n < n_cuts; ++n){
        cut_flags.push_back(true);
    }
    return cut_flags;
}

std::vector<bool> Check_Cuts(std::vector<double> cut_values, double lepton_property, std::vector<bool> cut_flags){
    for(int i = 0; i < cut_values.size(); ++i){
        if(lepton_property < cut_values[i]){
            cut_flags[i] = false;
        }
    }
    return cut_flags;
}

void Fill_Histogram(std::vector<TH1D*> h_varycuts, std::vector<bool> cut_flags, double reco_Higgs){
    for(int i = 0; i < cut_flags.size(); ++i){
        if(cut_flags[i] == true){
            h_varycuts[i] -> Fill(reco_Higgs);
        }
    }
}

void Write_Histogram(std::vector<TH1D*> h_varycuts){
    for(int i = 0; i < h_varycuts.size(); ++i){
        h_varycuts[i]->Write();
    }
}

