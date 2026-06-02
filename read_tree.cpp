#include "read_tree.hpp"

bool ReadTree::init(std::string filename)
{
    TFile* fFile = TFile::Open(filename.c_str(),"read");
    if (!fFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return true;
    }
    TTree* tree = (TTree*)fFile->Get("trackTree");
    if(!tree) {
        std::cerr << "trackTree not found in file\n";
        fFile->Close();
        return true;
    }

    fChain = tree;
    fChain->SetMakeClass(1);

    chg = nullptr;
    px = nullptr;
    py = nullptr;
    pz = nullptr;
    pid = nullptr;

    genJetEta = nullptr;
    genJetPt = nullptr;
    genJetPhi = nullptr;
    genJetChargedMultiplicity = nullptr;

    genDau_chg = nullptr;
    genDau_pid = nullptr;
    genDau_pt = nullptr;
    genDau_eta = nullptr;
    genDau_phi = nullptr;

    fChain->SetBranchAddress("chg", &chg);
    fChain->SetBranchAddress("px", &px);
    fChain->SetBranchAddress("py", &py);
    fChain->SetBranchAddress("pz", &pz);
    fChain->SetBranchAddress("pid", &pid);

    fChain->SetBranchAddress("genJetEta", &genJetEta);
    fChain->SetBranchAddress("genJetPt", &genJetPt);
    fChain->SetBranchAddress("genJetPhi", &genJetPhi);
    fChain->SetBranchAddress("genJetChargedMultiplicity", &genJetChargedMultiplicity);

    fChain->SetBranchAddress("genDau_chg", &genDau_chg);
    fChain->SetBranchAddress("genDau_pid", &genDau_pid);
    fChain->SetBranchAddress("genDau_pt", &genDau_pt);
    fChain->SetBranchAddress("genDau_eta", &genDau_eta);
    fChain->SetBranchAddress("genDau_phi", &genDau_phi);

    nEntries = fChain->GetEntriesFast();
    return false;
}
bool ReadTree::readEntry(Long64_t entry)
{
    if (!fChain) return false;
    fChain->LoadTree(entry);
    return fChain->GetEntry(entry) > 0;
}
void ReadTree::close()
{
    if (fChain) {
        TFile* fFile = fChain->GetCurrentFile();
        if (fFile) {
            fFile->Close();
        }
    }
}