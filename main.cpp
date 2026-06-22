#include "histograms.hpp"
#include "histograms.cpp"
#include "read_tree.cpp"      // add this line
#include "analysis.cpp"
#include <fstream>

auto readFileList = [](const char* path) {
    std::ifstream ifs(path);
    std::vector<std::string> files;
    std::string line;
    while(std::getline(ifs, line))
        if(!line.empty()) files.push_back(line);
    return files;
};

int main(int argc, char* argv[]) {
    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " input.txt input_wFSI.txt\n";
        return 1;
    }

    Histograms hist;
    fill_histograms(readFileList(argv[1]), hist);

    Histograms hist_fsi("FSI_");
    fill_histograms(readFileList(argv[2]), hist_fsi);

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

    hist.Draw(&hist_fsi);
    return 0;
}