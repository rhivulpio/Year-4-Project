#include <vector>

std::vector<double> Define_Cut_Values(int n_cuts, double min_cut, double step_cut);
std::vector<TString> Generate_Histogram_List(int n_pT_cuts, TString property);
void Make_Histograms(TString property, TString units, int n_cuts, double min_cut, double step_cut);
void Significance_Plots(TString property, TString units, int n_cuts, double min_cut, double step_cut);

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
        h_varycuts.push_back("Mass Reco Cuts/h_" + property + suffix + std::to_string(i));
    }
    return h_varycuts;
}