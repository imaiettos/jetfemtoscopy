#pragma once

#include "config.hpp"

class ReadTree {
public:
    ReadTree() = default;

    bool init(std::string filename);
    bool readEntry(Long64_t entry);
    void close();
    // Branch data
    std::vector<int>* chg = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<int>* pid = nullptr;

    std::vector<float>* genJetEta = nullptr;
    std::vector<float>* genJetPt = nullptr;
    std::vector<float>* genJetPhi = nullptr;
    std::vector<int>* genJetChargedMultiplicity = nullptr;

    std::vector<std::vector<int>>* genDau_chg = nullptr;
    std::vector<std::vector<int>>* genDau_pid = nullptr;
    std::vector<std::vector<float>>* genDau_pt = nullptr;
    std::vector<std::vector<float>>* genDau_eta = nullptr;
    std::vector<std::vector<float>>* genDau_phi = nullptr;

    Long64_t nEntries = 0;

private:
    TTree* fChain = nullptr;
};