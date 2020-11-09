
#include "Process.h"

bool Debug = true;

//defining the selection criteria for each particle
bool MuonCut(GenParticle * muon){
    double eta_min = -4.6; //eta range of the inner tracker
    double eta_max = 5.3;
    double pT_min = 15;
    if (eta_min<muon->Eta && muon->Eta<eta_max && muon->PT>pT_min){
        return true;
    }
    return false;
}

bool NeutrinoCut(GenParticle * neutrino){
    double eta_min = -1.9; //eta range of the EM calorimeter
    double eta_max = 2.4;
    double pT_min = 15;
    if (eta_min<neutrino->Eta && neutrino->Eta<eta_max && neutrino->PT>pT_min){
        return true;
    }
    return false;
}

bool JetCut(Jet * quarkjet){
    double eta_min = -5.0; //eta range of the hadronic calorimeter
    double eta_max = 5.5;
    double pT_min = 15;
    if (eta_min<quarkjet->Eta && quarkjet->Eta<eta_max && quarkjet->PT>pT_min){
        return true;
    }
    return false;
}

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

    h_EventCount = new TH1D("h_EventCount","",1,0,1);
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
    h_nu_eta = new TH1D("h_nu_eta","; Muon Neutrino Eta ; Events", 50, -6.0, 4.0);
    h_Higgs_eta = new TH1D("h_Higgs_eta","; Higgs Eta ; Events", 100, -8.0, 0); 

    //transverse energy histograms

    h_Jet_Et = new TH1D("h_Jet_Et","; Jet Transverse Energy [GeV]; Events", 100, 0, 160);
    h_Z_Et = new TH1D("h_Z_Et","; Z Transverse Energy [GeV]; Events", 100, 0, 180);
    h_mu_Et = new TH1D("h_mu_Et", "; Muon Transverse Energy [GeV] ; Events", 100, 0, 140);
    h_nu_Et = new TH1D("h_nu_Et","; Muon Neutrino Transverse Energy [GeV] ; Events", 100, 0, 200);
    h_Higgs_Et = new TH1D("h_Higgs_Et","; Higgs Transverse Energy [GeV]; Events", 100, 100, 240);

    //------------------------------------

    // Run the selection
    Process(reader);

    std::cout << "Events in EventCount: " << h_EventCount->GetEntries() << std::endl;

    std::cout << "Write to file..." << std::endl;

    OutputFile->cd();

    h_EventCount->Write();
    h_WeightCount->Write();

    h_Z_Pt->Write();
    h_Lepton_Pt->Write();
    h_Jet_Pt->Write();
    h_ZZ_Mass->Write();
    h_mu_pT->Write();

    h_mu_eta->Write();
    h_Jet_eta->Write();
    h_Z_eta->Write();
    h_nu_eta->Write();

    h_mu_Et->Write();
    h_Jet_Et->Write();
    h_Z_Et->Write();
    h_nu_Et->Write();

    h_Higgs_pT->Write();
    h_Higgs_eta->Write();
    h_Higgs_Et->Write();

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

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        HepMCEvent * event = (HepMCEvent*) bEvent->At(0);
    	const float Event_Weight = event->Weight;

        h_EventCount->Fill(0.5);
        h_WeightCount->Fill(0.5,Event_Weight);

        if( (entry > 0 && entry%1000 == 0) || Debug) {
            std::cout << "-------------------------------------------------------------"  << std::endl;
            std::cout << "Processing Event Number =  " << entry  << std::endl;
            std::cout << "-------------------------------------------------------------"  << std::endl;
        }

        //------------------------------------------------------------------
        // Jet Loop
        //------------------------------------------------------------------

        for(int i = 0; i < bJet->GetEntriesFast(); ++i) {

            Jet* jet = (Jet*) bJet->At(i);

            if(Debug) std::cout << "Jet " << i << " pT = " << jet->PT << " eta = " << jet->Eta << " phi = " << jet->Phi << 
            " mass = " << jet->Mass << " flavour = " << jet->Flavor << std::endl;

            TLorentzVector Vec_Jet;
            Vec_Jet.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);

            if(JetCut(jet)){
                h_Jet_Pt->Fill( Vec_Jet.Pt(), Event_Weight );
                h_Jet_eta->Fill(Vec_Jet.Eta(), Event_Weight);
                h_Jet_Et->Fill(TMath::Sqrt(Vec_Jet.Pt() * Vec_Jet.Pt() + Vec_Jet.M() * Vec_Jet.M()), Event_Weight);
            }

        } // Jet Loop


        //------------------------------------------------------------------
        // Lepton Loop
        //------------------------------------------------------------------
        4mu = true; // to check if the event satisfies the 4mu subchannel

        for(int i = 0; i < bTruthLepton->GetEntriesFast(); ++i) {

            GenParticle* lep = (GenParticle*) bTruthLepton->At(i);

            if(Debug) std::cout << "Lepton " << i << " PID = " << lep->PID << " pT = " << lep->PT << 
            " << eta = " << lep->Eta << " phi = " << lep->Phi << " mass = " << lep->Mass << 
            std::endl;

            
            TLorentzVector Vec_Lepton;
            Vec_Lepton.SetPtEtaPhiM(lep->PT,lep->Eta,lep->Phi,lep->Mass);
            

            // Look for electrons or muons
            if( abs(lep->PID) == 11 || abs(lep->PID) == 13 ) {

                    h_Lepton_Pt->Fill( Vec_Lepton.Pt(), Event_Weight );  

            }

            //looking for muons
            if( abs(lep->PID) == 13 && MuonCut(lep)) {
                h_mu_eta -> Fill(lep->Eta, Event_Weight);
                h_mu_Et -> Fill (TMath::Sqrt(lep->PT * lep->PT + lep->Mass * lep->Mass), Event_Weight);
                h_mu_pT-> Fill(Vec_Lepton.Pt(), Event_Weight);
            }

            //looking for neutrinos
            if( abs(lep->PID) == 12 && NeutrinoCut(lep) ){
                h_nu_Et -> Fill(TMath::Sqrt(lep->PT * lep->PT + lep->Mass * lep->Mass), Event_Weight);
                h_nu_eta -> Fill(lep->Eta, Event_Weight);
            }

            if(abs(lep->PID) == 11){ //if any electrons are seen in the decay then it is not a 4mu decay mode and so will be ignored
                4mu = false;
            }

        } // Lepton Loop

        //another lepton loop but only including the 4mu decays
        

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

            //looking for Z bosons
            if( abs(boson->PID) == 23 ){
		        h_Z_Pt->Fill( Vec_Boson.Pt(), Event_Weight);
                h_Z_eta->Fill(Vec_Boson.Eta(), Event_Weight);
                h_Z_Et->Fill(TMath::Sqrt(Vec_Boson.Pt() * Vec_Boson.Pt() + Vec_Boson.M() * Vec_Boson.M()), Event_Weight);
                list_Zboson.push_back(Vec_Boson);
            }

            //looking for Higgs bosons
            if( abs(boson->Mass) == 125 ){
		        h_Higgs_pT->Fill(Vec_Boson.Pt(), Event_Weight);
                h_Higgs_eta->Fill(Vec_Boson.Eta(), Event_Weight);
                h_Higgs_Et->Fill(TMath::Sqrt(Vec_Boson.Pt() * Vec_Boson.Pt() + Vec_Boson.M() * Vec_Boson.M()), Event_Weight);
        }
        }
        //Z/W/H Boson Loop

	if( list_Zboson.size() > 1 ) {

		TLorentzVector Higgs = list_Zboson.at(0) + list_Zboson.at(1);

              h_ZZ_Mass->Fill( Higgs.M() , Event_Weight );
	}


    } // Loop over all events
}