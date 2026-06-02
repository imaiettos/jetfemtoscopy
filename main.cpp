#include "histograms.hpp"
#include "histograms.cpp"
#include "read_tree.cpp"      // add this line
#include "analysis.cpp"
#include <fstream>

int main(int argc, char* argv[]) {
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " filelist.txt\n";
        return 1;
    }
    std::ifstream ifs(argv[1]);
    std::vector<std::string> files;
    std::string line;
    while(std::getline(ifs, line))
        if(!line.empty()) files.push_back(line);

    Histograms hist;
    fill_histograms(files, hist);

    TFile* out = TFile::Open("output_new.root", "RECREATE");
    for(int i=0; i<jet_nch_bin_count; i++)
        for(int j=0; j<jet_pt_bin_count; j++)
            for(int k=0; k<jet_eta_bin_count; k++){
                hist.hSignal_same[i][j][k]->Write();
                hist.hSignal_opposite[i][j][k]->Write();
                hist.hBckrnd_same[i][j][k]->Write();
                hist.hBckrnd_opposite[i][j][k]->Write();
                hist.hCorrelation_same[i][j][k]->Write();
                hist.hCorrelation_opposite[i][j][k]->Write();
                hist.hRelative[i][j][k]->Write();
            }
    out->Close();

    hist.Draw("correlations.pdf", "relative.pdf");
    return 0;
}