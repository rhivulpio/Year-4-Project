#include "TFile.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TVector.h"
#include "TError.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TVector3.h"

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"
#include "TRandom3.h"
#include "TClonesArray.h"

#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include "TEfficiency.h"

// Plots

TClonesArray * bEvent;
TClonesArray * bTruthLepton;
TClonesArray * bTruthWZ;
TClonesArray * bJet;


// Output
TFile * OutputFile;

TH1D * h_EventCount;
TH1D * h_WeightCount;

TH1D * h_Z_Pt;
TH1D * h_Lepton_Pt;
TH1D * h_Jet_Pt;
TH1D * h_ZZ_Mass;

//my graphs

//pseudorapidity and transverse momentum
TH1D * h_mu_pT;
TH1D * h_mu_eta;
TH1D * h_Jet_eta; 
TH1D * h_Z_eta; 
TH1D * h_nu_eta; 

//transverse energy
TH1D * h_mu_Et;
TH1D * h_Jet_Et;
TH1D * h_Z_Et; 
TH1D * h_nu_Et;

//Higgs graphs
TH1D * h_Higgs_pT;
TH1D * h_Higgs_eta;
TH1D * h_Higgs_Et;

//4mu event graphs
TH1D * h4mu_mu_pT;
TH1D * h4mu_mu_Et;
TH1D * h4mu_nu_Et;
TH1D * h4mu_mu_eta;
TH1D * h4mu_nu_eta; 

//acceptance plots
TEfficiency * e_eta;
TEfficiency * e_Et;
TEfficiency * e_pT;
TEfficiency * e4mu_eta;
TEfficiency * e4mu_Et;
TEfficiency * e4mu_pT;

TEfficiency * e_H_eta;
TEfficiency * e_H_Et;
TEfficiency * e_H_pT;

TEfficiency * e4mu_H_eta;
TEfficiency * e4mu_H_Et;
TEfficiency * e4mu_H_pT;

//plotting variables against each other
TEfficiency * e_mu_pT_eta;

//truejet histograms
TH1D * h_trueJet_Et;
TH1D * h_trueJet_eta;
TH1D * h_trueJet_Pt;

//2D missing energy plot
TH2D * h_ME_nu_Et;
TH2D * h_ME_nu_eta;
TH2D * h_ME_nu_phi;

TH2D * h_mu_pT_eta;

//kinematic reconstruction
TH2D * x_Qsquared_electron;
TH2D * x_Qsquared_hadron;
TH2D * log_Qsquared_plot;
TH2D * log_x_plot;
TH2D * log_y_plot;
TH1D * h_logx_electron;
TH1D * h_logy_electron;
TH1D * h_logQsquared_electron;
TH1D * h_logx_hadron;
TH1D * h_logy_hadron;
TH1D * h_logQsquared_hadron;

ExRootTreeReader * InitReader(const TString FilePath);

void Process(ExRootTreeReader * treeReader);

void ClearBranches();

int main(int argc, char* argv[]);
