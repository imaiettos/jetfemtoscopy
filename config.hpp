#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
using TMath::ATan;
using TMath::Exp;

constexpr float pion_mass = 0.13957f;

constexpr int jet_nch_bin_count = 10;
constexpr int jet_pt_bin_count = 4;
constexpr int jet_eta_bin_count = 1;

constexpr float jet_pt_cut  = 100.0f;
constexpr float jet_eta_cut = 1.6f;

constexpr float Q_min       = 0.0f;   // start from 0; underflow handled by Q_lower_cutoff
//constexpr float Q_max       = 2.0f;   // 2 GeV is more than enough for HBT
//constexpr int   Q_bin_count = 200;    // now each bin is 10 MeV — much better resolution

constexpr int Q_lower_cutoff = 2;

constexpr float Q_max       = 2.0f;
constexpr int   Q_bin_count = 100;  // keep 10 MeV bins

//constexpr float Q_norm_min = 4.0f;
//constexpr float Q_norm_max = 6.0f;

constexpr float Q_norm_min = 0.8f;
constexpr float Q_norm_max = 1.5f;


inline const double jet_nch_low[jet_nch_bin_count] = {0,0,25,36,48,60,71,78,91,97};
inline const double jet_nch_high[jet_nch_bin_count] = {1000,25,36,48,60,71,78,91,97,1000};

inline const double jet_pt_low[jet_pt_bin_count]  = {500, 500, 620, 800};
inline const double jet_pt_high[jet_pt_bin_count] = {2000, 620, 800, 2000};

inline const double jet_eta_low[jet_eta_bin_count] = {-1.6};
inline const double jet_eta_high[jet_eta_bin_count] = {1.6};
