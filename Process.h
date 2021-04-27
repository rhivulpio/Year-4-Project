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
#include "THStack.h"

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
TH1D * h4mu_nu_pT;
TH1D * h4mu_mu_eta;
TH1D * h4mu_nu_eta; 

//4mu event seen by detector graphs
TH1D * h4mu_mu_pT_seen;
TH1D * h4mu_mu_Et_seen;
TH1D * h4mu_nu_pT_seen;
TH1D * h4mu_mu_eta_seen;
TH1D * h4mu_nu_eta_seen; 
TH1D * h4mu_Z_pT_seen;
TH1D * h4mu_Z_eta_seen;
TH1D * h4mu_Higgs_pT_seen;
TH1D * h4mu_Higgs_eta_seen;
TH1D * h4mu_jet_eta_seen;
TH1D * h4mu_jet_pT_seen;

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
TH2D * h_ME_nu_pT;
TH2D * h_ME_nu_eta;
TH2D * h_ME_nu_phi;

//scatter plots
TH2D * h_mu_pT_eta;
TH2D * h_4mu_pT_eta;

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

//Z mass reconstruction
TH1D * h_Higgs_reco;
TH1D * h_ZZ_mass_reco;
TH1D * h_Z_reco;
TH1D * h_Zstar_reco;

//event weight
TH1D * h_eventweight;

//mass reconstruction cuts
TH1D * h_new;
std::vector<TH1D*> h_pT_cuts;
std::vector<TH1D*> h_Zstar_cuts;
std::vector<TH1D*> h_Z_cuts;

ExRootTreeReader * InitReader(const TString FilePath);

void Process(ExRootTreeReader * treeReader); //removed bool signal from the arguments

void ClearBranches();

int main(int argc, char* argv[]);

std::vector<double> Electron_Reconstruction(TLorentzVector nu);
std::vector<double> Hadron_Reconstruction(TLorentzVector Missing_Energy_Vector);
std::vector<TH1D*> Define_Histograms(TString hist_name, int n_cuts);
std::vector<bool> Initialise_Flags(int n_cuts);
std::vector<bool> Check_Cuts(std::vector<double> cut_values, double lepton_property, std::vector<bool> cut_flags);
void Fill_Histogram(std::vector<TH1D*> h_varycuts, std::vector<bool> cut_flags, double reco_Higgs);
void Write_Histogram(std::vector<TH1D*> h_varycuts);
std::vector<TLorentzVector> Lorentz_Vector(std::vector<GenParticle*> particles);
void Particles_Antiparticles(std::vector<GenParticle*> all_muons_seen, std::vector<GenParticle*> all_electrons_seen, std::vector<GenParticle*> &particles, std::vector<GenParticle*> &antiparticles);
std::vector<double> Mass_Reconstruction(std::vector<GenParticle*> all_muons_seen, std::vector<GenParticle*> all_electrons_seen, std::vector<TLorentzVector> particles, std::vector<TLorentzVector> antiparticles);
std::tuple<std::vector<TLorentzVector>, std::vector<int>> Smear(std::vector<GenParticle*> particles);
