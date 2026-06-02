#include "read_tree.hpp"
#include "histograms.hpp"
#include "cuts.cpp"
#include "physics.cpp"

struct ParticleInfo {
    ROOT::Math::PtEtaPhiMVector p4;
    ROOT::Math::PtEtaPhiMVector jet4v;
    bool charge;

    ParticleInfo(ROOT::Math::PtEtaPhiMVector _p4, ROOT::Math::PtEtaPhiMVector _jet4v, int _charge)
        : p4(_p4), jet4v(_jet4v), charge(_charge) {}
};

void fill_histograms(const std::vector<std::string>& fileList, Histograms& hist){
    ReadTree reader;
    std::deque<std::vector<ParticleInfo>> jetPool[jet_nch_bin_count][jet_pt_bin_count][jet_eta_bin_count];
    unsigned int max_pool_size = 50;

    for(int f = 0; f<(int)fileList.size(); f++){
        if(reader.init(fileList.at(f).c_str())) {
            std::cerr << "Error initializing reader for file: " << fileList.at(f) << std::endl;
            continue;
        }

        for (Long64_t ievent=0; ievent < reader.nEntries; ievent++){
            reader.readEntry(ievent);
            if(!eventPass(reader.genJetPt, reader.genJetPt->size())) continue;

            for(size_t ijet=0; ijet < reader.genJetPt->size(); ijet++){
                if(!jetPass(reader.genJetEta, reader.genJetPt, ijet)) continue;

                float jet_pt  = reader.genJetPt->at(ijet);
                float jet_eta = reader.genJetEta->at(ijet);
                int   jet_nch = reader.genJetChargedMultiplicity->at(ijet);

                std::vector<int> nch_bins;
                for(int i = 0; i < jet_nch_bin_count; i++)
                    if(jet_nch >= jet_nch_low[i] && jet_nch < jet_nch_high[i])
                        nch_bins.push_back(i);

                std::vector<int> pt_bins;
                for(int j = 0; j < jet_pt_bin_count; j++)
                    if(jet_pt >= jet_pt_low[j] && jet_pt < jet_pt_high[j])
                        pt_bins.push_back(j);

                std::vector<int> eta_bins;
                for(int k = 0; k < jet_eta_bin_count; k++)
                    if(jet_eta >= jet_eta_low[k] && jet_eta < jet_eta_high[k])
                        eta_bins.push_back(k);

                size_t num_dau = (reader.genDau_pt->at(ijet)).size();

                // reconstruct jet 4-vector from daughters
                ROOT::Math::PtEtaPhiMVector jet4v(0,0,0,0);
                for(size_t i=0; i<num_dau; i++){
                    ROOT::Math::PtEtaPhiMVector dau(
                        (reader.genDau_pt->at(ijet)).at(i),
                        (reader.genDau_eta->at(ijet)).at(i),
                        (reader.genDau_phi->at(ijet)).at(i),
                        0.13957f);
                    jet4v = jet4v + dau;
                }

                std::vector<ParticleInfo> current_jet;

                for(size_t i=0; i<num_dau; i++){
                    if(std::abs((reader.genDau_pid->at(ijet)).at(i)) != 211) continue;
                    if(!nchPass((reader.genDau_eta->at(ijet)).at(i), (reader.genDau_pt->at(ijet)).at(i))) continue;

                    ROOT::Math::PtEtaPhiMVector p1(
                        (reader.genDau_pt->at(ijet)).at(i),
                        (reader.genDau_eta->at(ijet)).at(i),
                        (reader.genDau_phi->at(ijet)).at(i),
                        pion_mass);
                    current_jet.emplace_back(p1, jet4v, (reader.genDau_pid->at(ijet)).at(i) > 0);

                    // signal pairs — same jet
                    for(size_t j=i+1; j<num_dau; j++){
                        if(std::abs((reader.genDau_pid->at(ijet)).at(j)) != 211) continue;
                        if(!nchPass((reader.genDau_eta->at(ijet)).at(j), (reader.genDau_pt->at(ijet)).at(j))) continue;
                        ROOT::Math::PtEtaPhiMVector p2(
                            (reader.genDau_pt->at(ijet)).at(j),
                            (reader.genDau_eta->at(ijet)).at(j),
                            (reader.genDau_phi->at(ijet)).at(j),
                            pion_mass);
                        if(!pairPass(p1.Phi(), p1.Eta(), p2.Phi(), p2.Eta())) continue;
                        // BEFORE
                        //float Q_inv = GetQ_JetFrame(p1, p2, jet4v);

                        // AFTER
                        float Q_inv = GetQ_JetAxis(p1, p2, jet4v);
                        for(int nch_bin_num : nch_bins)
                            for(int pt_bin_num : pt_bins)
                                for(int eta_bin_num : eta_bins){
                                    if((reader.genDau_pid->at(ijet)).at(i) == (reader.genDau_pid->at(ijet)).at(j))
                                        hist.hSignal_same[nch_bin_num][pt_bin_num][eta_bin_num]->Fill(Q_inv);
                                    else
                                        hist.hSignal_opposite[nch_bin_num][pt_bin_num][eta_bin_num]->Fill(Q_inv);
                                }
                    }

                    // background pairs — mixed jets
                    for(int nch_bin_num : nch_bins){
                        for(int pt_bin_num : pt_bins){
                            for(int eta_bin_num : eta_bins){
                                for(const auto& past_jet : jetPool[nch_bin_num][pt_bin_num][eta_bin_num]){
                                    for(const auto& past_p : past_jet){
                                        if(!pairPass(p1.Phi(), p1.Eta(), past_p.p4.Phi(), past_p.p4.Eta())) continue;
                                        // BEFORE
                                        //float Q_inv = GetQ_JetFrame_Double(p1, jet4v, past_p.p4, past_p.jet4v);

                                        // AFTER
                                        float Q_inv = GetQ_JetAxis_Double(p1, jet4v, past_p.p4, past_p.jet4v);                                        if(past_p.charge == current_jet.back().charge)
                                            hist.hBckrnd_same[nch_bin_num][pt_bin_num][eta_bin_num]->Fill(Q_inv);
                                        else
                                            hist.hBckrnd_opposite[nch_bin_num][pt_bin_num][eta_bin_num]->Fill(Q_inv);
                                    }
                                }
                            }
                        }
                    }
                } // end particle loop

                if(current_jet.size() > 0){
                    for(int nch_bin_num : nch_bins)
                        for(int pt_bin_num : pt_bins)
                            for(int eta_bin_num : eta_bins){
                                jetPool[nch_bin_num][pt_bin_num][eta_bin_num].push_back(current_jet);
                                if(jetPool[nch_bin_num][pt_bin_num][eta_bin_num].size() > max_pool_size)
                                    jetPool[nch_bin_num][pt_bin_num][eta_bin_num].pop_front();
                            }
                }
            } // end jet loop
        } // end event loop
        reader.close();
    } // end file loop

    // calculating correlation functions
    for(int i = 0; i < jet_nch_bin_count; i++){
        float nch_bin_min = jet_nch_low[i];
        float nch_bin_max = jet_nch_high[i];
        for(int j = 0; j < jet_pt_bin_count; j++){
            float pt_bin_min = jet_pt_low[j];
            float pt_bin_max = jet_pt_high[j];
            for(int k = 0; k < jet_eta_bin_count; k++){
                float eta_bin_min = jet_eta_low[k];
                float eta_bin_max = jet_eta_high[k];
                TH1D* hS1 = hist.hSignal_same[i][j][k];
                TH1D* hB1 = hist.hBckrnd_same[i][j][k];
                TH1D* hS2 = hist.hSignal_opposite[i][j][k];
                TH1D* hB2 = hist.hBckrnd_opposite[i][j][k];
                double S1_integral = hS1->Integral(hist.Q_norm_min_bin, hist.Q_norm_max_bin);
                double B1_integral = hB1->Integral(hist.Q_norm_min_bin, hist.Q_norm_max_bin);
                double S2_integral = hS2->Integral(hist.Q_norm_min_bin, hist.Q_norm_max_bin);
                double B2_integral = hB2->Integral(hist.Q_norm_min_bin, hist.Q_norm_max_bin);
                if(B1_integral > 0) hB1->Scale(S1_integral/B1_integral);
                else hB1->Reset();
                if(B2_integral > 0) hB2->Scale(S2_integral/B2_integral);
                else hB2->Reset();
                TString hname_corr_same     = Form("Corr_Same_%d_%d_%d", i, j, k);
                TString hname_corr_opposite = Form("Corr_Opposite_%d_%d_%d", i, j, k);
                TString htitle_corr_same = Form(
                    "Correlation: Same Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    nch_bin_min, nch_bin_max, pt_bin_min, pt_bin_max, eta_bin_min, eta_bin_max);
                TString htitle_corr_opposite = Form(
                    "Correlation: Opposite Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    nch_bin_min, nch_bin_max, pt_bin_min, pt_bin_max, eta_bin_min, eta_bin_max);
                hist.hCorrelation_same[i][j][k] = (TH1D*)hS1->Clone(hname_corr_same);
                if(hB1->GetEntries()) hist.hCorrelation_same[i][j][k]->Divide(hS1, hB1, 1.0, 1.0, "B");
                else hist.hCorrelation_same[i][j][k]->Reset();
                hist.hCorrelation_same[i][j][k]->SetTitle(htitle_corr_same);
                hist.hCorrelation_opposite[i][j][k] = (TH1D*)hS2->Clone(hname_corr_opposite);
                if(hB2->GetEntries()) hist.hCorrelation_opposite[i][j][k]->Divide(hS2, hB2, 1.0, 1.0, "B");
                else hist.hCorrelation_opposite[i][j][k]->Reset();
                hist.hCorrelation_opposite[i][j][k]->SetTitle(htitle_corr_opposite);
                TString hname_relative  = Form("Relative_%d_%d_%d", i, j, k);
                TString htitle_relative = Form(
                    "Relative Correlation, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    nch_bin_min, nch_bin_max, pt_bin_min, pt_bin_max, eta_bin_min, eta_bin_max);
                hist.hRelative[i][j][k] = (TH1D*)hist.hCorrelation_same[i][j][k]->Clone(hname_relative);
                if(hist.hCorrelation_opposite[i][j][k]->GetEntries())
                    hist.hRelative[i][j][k]->Divide(hist.hCorrelation_opposite[i][j][k]);
                else hist.hRelative[i][j][k]->Reset();
                hist.hRelative[i][j][k]->SetTitle(htitle_relative);
                for(int q_bin = 0; q_bin < Q_lower_cutoff; q_bin++){
                    hist.hCorrelation_same[i][j][k]->SetBinContent(q_bin, 0.0);
                    hist.hCorrelation_same[i][j][k]->SetBinError(q_bin, 0.0);
                    hist.hCorrelation_opposite[i][j][k]->SetBinContent(q_bin, 0.0);
                    hist.hCorrelation_opposite[i][j][k]->SetBinError(q_bin, 0.0);
                    hist.hRelative[i][j][k]->SetBinContent(q_bin, 0.0);
                    hist.hRelative[i][j][k]->SetBinError(q_bin, 0.0);
                }
            }
        }
    }
}