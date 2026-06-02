#pragma once

#include "config.hpp"

class Histograms {
public:
    Histograms();
    void Draw(std::string filename1, std::string filename2);
    
    TH1D* hSignal_same[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hSignal_opposite[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hBckrnd_same[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hBckrnd_opposite[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hCorrelation_same[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hCorrelation_opposite[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    TH1D* hRelative[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    int Q_norm_min_bin;
    int Q_norm_max_bin;
};