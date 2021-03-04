#include <vector>

std::vector<double> Define_Cut_Values(int n_cuts, double min_cut, double step_cut);

std::vector<double> Define_Cut_Values(int n_cuts, double min_cut, double step_cut){
    std::vector<double> cut_values;
	for(int n = 0; n < n_cuts; ++n){	
        cut_values.push_back(min_cut + n*step_cut);
	}
    return cut_values;
}