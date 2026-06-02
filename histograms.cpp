#include "histograms.hpp"
#include "config.hpp"

Histograms::Histograms()
{
    TH1::SetDefaultSumw2(kTRUE);

    // NEW: zero out derived histogram pointers before fill_histograms sets them
    for(int i=0; i<jet_nch_bin_count; i++)
        for(int j=0; j<jet_pt_bin_count; j++)
            for(int k=0; k<jet_eta_bin_count; k++){
                hCorrelation_same[i][j][k]     = nullptr;
                hCorrelation_opposite[i][j][k] = nullptr;
                hRelative[i][j][k]             = nullptr;
            }

    for(int i=0; i<jet_nch_bin_count; i++){
        for(int j=0; j<jet_pt_bin_count; j++){
            for(int k=0; k<jet_eta_bin_count; k++){
                TString hname_signal_same = Form("Signal_Same_%d_%d_%d", i, j, k);
                TString hname_signal_opposite = Form("Signal_Opposite_%d_%d_%d", i, j, k);
                TString hname_bkg_same = Form("Background_Same_%d_%d_%d", i, j, k);
                TString hname_bkg_opposite = Form("Background_Opposite_%d_%d_%d", i, j, k);
                TString htitle_signal_same = Form(
                    "Signal: Same Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    jet_nch_low[i], jet_nch_high[i],
                    jet_pt_low[j],  jet_pt_high[j],
                    jet_eta_low[k], jet_eta_high[k]
                );
                TString htitle_signal_opposite = Form(
                    "Signal: Opposite Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    jet_nch_low[i], jet_nch_high[i],
                    jet_pt_low[j],  jet_pt_high[j],
                    jet_eta_low[k], jet_eta_high[k]
                );
                TString htitle_bkg_same = Form(
                    "Background: Same Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    jet_nch_low[i], jet_nch_high[i],
                    jet_pt_low[j],  jet_pt_high[j],
                    jet_eta_low[k], jet_eta_high[k]
                );
                TString htitle_bkg_opposite = Form(
                    "Background: Opposite Charge, %.0f < N_{ch} < %.0f, %.1f < p_{T} < %.1f, %.2f < #eta < %.2f",
                    jet_nch_low[i], jet_nch_high[i],
                    jet_pt_low[j],  jet_pt_high[j],
                    jet_eta_low[k], jet_eta_high[k]
                );
                hSignal_same[i][j][k] = new TH1D(hname_signal_same, htitle_signal_same, Q_bin_count, Q_min, Q_max);
                hSignal_opposite[i][j][k] = new TH1D(hname_signal_opposite, htitle_signal_opposite, Q_bin_count, Q_min, Q_max);
                hBckrnd_opposite[i][j][k] = new TH1D(hname_bkg_opposite, htitle_bkg_opposite, Q_bin_count, Q_min, Q_max);
                hBckrnd_same[i][j][k] = new TH1D(hname_bkg_same, htitle_bkg_same, Q_bin_count, Q_min, Q_max);
            }
        }
    }
    Q_norm_min_bin = hSignal_same[0][0][0]->FindBin(Q_norm_min);
    Q_norm_max_bin = hSignal_same[0][0][0]->FindBin(Q_norm_max);
}

void Histograms::Draw(std::string filename1, std::string filename2)
{
    for(int i=0; i<jet_nch_bin_count; i++)
        for(int j=0; j<jet_pt_bin_count; j++)
            for(int k=0; k<jet_eta_bin_count; k++)
                if(!hCorrelation_same[i][j][k] || !hCorrelation_opposite[i][j][k] || !hRelative[i][j][k]){
                    std::cerr << "Draw() called before fill_histograms completed.\n";
                    return;
                }

    TCanvas* c = new TCanvas("c", "Histograms", 800, 600);
    TLine line(Q_min, 1.0, Q_max, 1.0);
    line.SetLineColor(kBlack);
    line.SetLineStyle(7);
    line.SetLineWidth(1);

    c->Print((filename1 + "[").c_str());  // open PDF
    for(int i = 0; i < jet_nch_bin_count; i++){
        for(int j = 0; j < jet_pt_bin_count; j++){
            for(int k = 0; k < jet_eta_bin_count; k++){
                c->Clear();
                hCorrelation_same[i][j][k]->SetMarkerStyle(20);
                hCorrelation_same[i][j][k]->SetMarkerColor(kRed);
                hCorrelation_same[i][j][k]->SetLineColor(kRed);
                hCorrelation_same[i][j][k]->SetFillStyle(0);
                hCorrelation_same[i][j][k]->GetXaxis()->SetTitle("Q");
                hCorrelation_same[i][j][k]->GetYaxis()->SetTitle("C(Q)");
                hCorrelation_same[i][j][k]->GetXaxis()->SetRangeUser(Q_min, Q_max);
                hCorrelation_same[i][j][k]->Draw("E1");
                hCorrelation_opposite[i][j][k]->SetMarkerStyle(21);
                hCorrelation_opposite[i][j][k]->SetMarkerColor(kBlue);
                hCorrelation_opposite[i][j][k]->SetLineColor(kBlue);
                hCorrelation_opposite[i][j][k]->SetFillStyle(0);
                hCorrelation_opposite[i][j][k]->Draw("E1 SAME");
                line.Draw("SAME");
                TLegend leg(0.65, 0.7, 0.88, 0.88);
                leg.AddEntry(hCorrelation_same[i][j][k], "Same charge", "p");
                leg.AddEntry(hCorrelation_opposite[i][j][k], "Opposite charge", "p");
                leg.Draw();
                c->Print(filename1.c_str());  // add page
            }
        }
    }
    c->Print((filename1 + "]").c_str());  // close PDF

    c->Print((filename2 + "[").c_str());  // open PDF
    for(int i = 0; i < jet_nch_bin_count; i++){
        for(int j = 0; j < jet_pt_bin_count; j++){
            for(int k = 0; k < jet_eta_bin_count; k++){
                c->Clear();
                hRelative[i][j][k]->SetMarkerStyle(20);
                hRelative[i][j][k]->SetMarkerColor(kGreen+2);
                hRelative[i][j][k]->SetLineColor(kGreen+2);
                hRelative[i][j][k]->SetFillStyle(0);
                hRelative[i][j][k]->GetXaxis()->SetTitle("Q");
                hRelative[i][j][k]->GetYaxis()->SetTitle("C_{same}(Q) / C_{opposite}(Q)");
                hRelative[i][j][k]->GetXaxis()->SetRangeUser(Q_min, Q_max);
                hRelative[i][j][k]->Draw("E1");
                line.Draw("SAME");
                c->Print(filename2.c_str());  // add page
            }
        }
    }
    c->Print((filename2 + "]").c_str());  // close PDF

    delete c;
}