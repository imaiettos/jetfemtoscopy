#include "histograms.hpp"
#include "config.hpp"
#include <sys/stat.h>

Histograms::Histograms(const std::string& prefix)
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
                TString hname_signal_same     = Form("%sSignal_Same_%d_%d_%d",       prefix.c_str(), i, j, k);
                TString hname_signal_opposite = Form("%sSignal_Opposite_%d_%d_%d",   prefix.c_str(), i, j, k);
                TString hname_bkg_same        = Form("%sBackground_Same_%d_%d_%d",   prefix.c_str(), i, j, k);
                TString hname_bkg_opposite    = Form("%sBackground_Opposite_%d_%d_%d", prefix.c_str(), i, j, k);
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

void Histograms::Draw(Histograms* fsi)
{
    for(int i=0; i<jet_nch_bin_count; i++)
        for(int j=0; j<jet_pt_bin_count; j++)
            for(int k=0; k<jet_eta_bin_count; k++)
                if(!hCorrelation_same[i][j][k] || !hCorrelation_opposite[i][j][k] || !hRelative[i][j][k]){
                    std::cerr << "Draw() called before fill_histograms completed.\n";
                    return;
                }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas* c = new TCanvas("c", "Histograms", 700, 700);
    c->SetMargin(0.15, 0.05, 0.15, 0.05);

    TLine line(Q_min, 1.0, Q_max, 1.0);
    line.SetLineColor(kBlack);
    line.SetLineStyle(7);
    line.SetLineWidth(1);

    mkdir("output", 0755);

    for(int i = 0; i < jet_nch_bin_count; i++){
        for(int j = 0; j < jet_pt_bin_count; j++){
            for(int k = 0; k < jet_eta_bin_count; k++){
                if(hSignal_same[i][j][k]->GetEntries() == 0 &&
                   hSignal_opposite[i][j][k]->GetEntries() == 0)
                    continue;

                TLatex tex;
                tex.SetNDC();
                tex.SetTextSize(0.035);
                tex.SetTextFont(42);

                // --- correlations ---
                c->Clear();
                hCorrelation_same[i][j][k]->SetMarkerStyle(20);
                hCorrelation_same[i][j][k]->SetMarkerColor(kBlack);
                hCorrelation_same[i][j][k]->SetLineColor(kBlack);
                hCorrelation_same[i][j][k]->SetFillStyle(0);
                hCorrelation_same[i][j][k]->GetXaxis()->SetTitle("Q");
                hCorrelation_same[i][j][k]->GetYaxis()->SetTitle("C(Q)");
                hCorrelation_same[i][j][k]->GetXaxis()->SetRangeUser(Q_min, Q_max);
                hCorrelation_same[i][j][k]->GetYaxis()->SetRangeUser(0.75, 2.05);
                hCorrelation_same[i][j][k]->Draw("E1");
                hCorrelation_opposite[i][j][k]->SetMarkerStyle(24);
                hCorrelation_opposite[i][j][k]->SetMarkerColor(kBlue);
                hCorrelation_opposite[i][j][k]->SetLineColor(kBlue);
                hCorrelation_opposite[i][j][k]->SetFillStyle(0);
                hCorrelation_opposite[i][j][k]->Draw("E1 SAME");
                if(fsi){
                    fsi->hCorrelation_same[i][j][k]->SetMarkerStyle(20);
                    fsi->hCorrelation_same[i][j][k]->SetMarkerColor(kRed);
                    fsi->hCorrelation_same[i][j][k]->SetLineColor(kRed);
                    fsi->hCorrelation_same[i][j][k]->SetFillStyle(0);
                    fsi->hCorrelation_same[i][j][k]->Draw("E1 SAME");
                    fsi->hCorrelation_opposite[i][j][k]->SetMarkerStyle(24);
                    fsi->hCorrelation_opposite[i][j][k]->SetMarkerColor(kOrange+1);
                    fsi->hCorrelation_opposite[i][j][k]->SetLineColor(kOrange+1);
                    fsi->hCorrelation_opposite[i][j][k]->SetFillStyle(0);
                    fsi->hCorrelation_opposite[i][j][k]->Draw("E1 SAME");
                }
                line.Draw("SAME");
                TLegend leg(0.5, 0.62, 0.88, 0.88);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.032);
                leg.AddEntry(hCorrelation_same[i][j][k],     "No FSI (SS)", "lp");
                leg.AddEntry(hCorrelation_opposite[i][j][k], "No FSI (OS)", "lp");
                if(fsi){
                    leg.AddEntry(fsi->hCorrelation_same[i][j][k],     "FSI (SS)", "lp");
                    leg.AddEntry(fsi->hCorrelation_opposite[i][j][k], "FSI (OS)", "lp");
                }
                leg.AddEntry(&line, "Reference (C=1)", "l");
                leg.Draw();
                tex.DrawLatex(0.18, 0.9,  Form("%.0f #leq N_{ch} < %.0f", jet_nch_low[i], jet_nch_high[i]));
                tex.DrawLatex(0.18, 0.85, Form("%.1f #leq p_{T} < %.1f GeV", jet_pt_low[j], jet_pt_high[j]));
                tex.DrawLatex(0.18, 0.8,  Form("%.2f #leq #eta < %.2f", jet_eta_low[k], jet_eta_high[k]));
                c->Print(Form("output/correlations_nch%d_pt%d_eta%d.pdf", i, j, k));

                // --- relative ---
                c->Clear();
                hRelative[i][j][k]->SetMarkerStyle(20);
                hRelative[i][j][k]->SetMarkerColor(kBlack);
                hRelative[i][j][k]->SetLineColor(kBlack);
                hRelative[i][j][k]->SetFillStyle(0);
                hRelative[i][j][k]->GetXaxis()->SetTitle("Q");
                hRelative[i][j][k]->GetYaxis()->SetTitle("C_{same}(Q) / C_{opposite}(Q)");
                hRelative[i][j][k]->GetXaxis()->SetRangeUser(Q_min, Q_max);
                hRelative[i][j][k]->GetYaxis()->SetRangeUser(0.75, 2.05);
                hRelative[i][j][k]->Draw("E1");
                TLegend leg2(0.5, 0.68, 0.88, 0.88);
                leg2.SetBorderSize(0);
                leg2.SetFillStyle(0);
                leg2.SetTextFont(42);
                leg2.SetTextSize(0.032);
                leg2.AddEntry(hRelative[i][j][k], "No FSI", "lp");
                if(fsi){
                    fsi->hRelative[i][j][k]->SetMarkerStyle(20);
                    fsi->hRelative[i][j][k]->SetMarkerColor(kRed);
                    fsi->hRelative[i][j][k]->SetLineColor(kRed);
                    fsi->hRelative[i][j][k]->SetFillStyle(0);
                    fsi->hRelative[i][j][k]->Draw("E1 SAME");
                    leg2.AddEntry(fsi->hRelative[i][j][k], "With FSI", "lp");
                }
                leg2.AddEntry(&line, "Reference (C=1)", "l");
                leg2.Draw();
                line.Draw("SAME");
                tex.DrawLatex(0.18, 0.9,  Form("%.0f #leq N_{ch} < %.0f", jet_nch_low[i], jet_nch_high[i]));
                tex.DrawLatex(0.18, 0.85, Form("%.1f #leq p_{T} < %.1f GeV", jet_pt_low[j], jet_pt_high[j]));
                tex.DrawLatex(0.18, 0.8,  Form("%.2f #leq #eta < %.2f", jet_eta_low[k], jet_eta_high[k]));
                c->Print(Form("output/relative_nch%d_pt%d_eta%d.pdf", i, j, k));
            }
        }
    }

    delete c;
}