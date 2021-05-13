#include <vector>

//------------------------------------
// Parameters for Selection Criteria
//------------------------------------
double pT_min = 5; //in GeV, used for the jets and leptons
double eta_min_j = -5.0; //eta range of the hadronic calorimeter, used for the jets
double eta_max_j = 5.6;
double eta_min_mu = -4.6; //eta range of the inner tracker, used for the muons
double eta_max_mu = 5.3;
double eta_min_e = -4.8; //eta range of the electromagnetic calorimeter, used for the electrons
double eta_max_e = 5.5;

//------------------------------------
// Parameters for Cuts Analysis
//------------------------------------
int n_pT_cuts = 40;
double min_pT_cut = 10.0; //GeV
double step_pT_cut = 2.0; //GeV

int n_Zstar_cuts = 40;
double min_Zstar_cut = 10.0; //GeV
double step_Zstar_cut = 1.0; //GeV

int n_Z_cuts = 40;
double min_Z_cut = 40.0; //GeV
double step_Z_cut = 2.0; //GeV

//------------------------------------
// Declaring Functions
//------------------------------------
std::vector<double> Define_Cut_Values(int n_cuts, double min_cut, double step_cut);
std::vector<TString> Generate_Histogram_List(int n_pT_cuts, TString property);
void Make_Histograms(TString property, TString units, int n_cuts, double min_cut, double step_cut);
void Significance_Plots(TString property, TString units, int n_cuts, double min_cut, double step_cut);
TCanvas * hstack(TString histogram);

std::vector<double> Define_Cut_Values(int n_cuts, double min_cut, double step_cut){
    std::vector<double> cut_values;
	for(int n = 0; n < n_cuts; ++n){	
        cut_values.push_back(min_cut + n*step_cut);
	}
    return cut_values;
}

std::vector<TString> Generate_Histogram_List(int n_cuts, TString property){
    std::vector<TString> h_varycuts;
    TString suffix = "_cut";
    for(int i = 0; i < n_cuts; ++i){
        h_varycuts.push_back("Vary Cuts/" + property + " Cuts/h_" + property + suffix + std::to_string(i));
    }
    return h_varycuts;
}