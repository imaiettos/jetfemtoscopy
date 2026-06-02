#pragma once

#include "config.hpp"


//WRONG
/*bool pairPass(float phi1, float eta1, float phi2, float eta2)
{
    return true;
    const float deltaEtaMin = 0.014;
    const float deltaPhiMin = 0.022;
    float deltaEta = eta1 - eta2;
    float deltaPhi = TMath::ACos(TMath::Cos(phi1 - phi2));
    if ((deltaEta/deltaEtaMin)*(deltaEta/deltaEtaMin) + (deltaPhi/deltaPhiMin)*(deltaPhi/deltaPhiMin) <= 1.0)
        return false;
    return true;
}
*/

// CORRECT
bool pairPass(float phi1, float eta1, float phi2, float eta2)
{
    const float deltaEtaMin = 0.014;
    const float deltaPhiMin = 0.022;
    float deltaEta = eta1 - eta2;
    float deltaPhi = TMath::ACos(TMath::Cos(phi1 - phi2));
    if ((deltaEta/deltaEtaMin)*(deltaEta/deltaEtaMin)
      + (deltaPhi/deltaPhiMin)*(deltaPhi/deltaPhiMin) <= 1.0)
        return false;
    return true;
}

bool nchPass(float eta, float pT)
{
    return (fabs(eta) < 2.4 && pT > 0.3);
}

//bool eventPass(std::vector< float > *jetPt, int jetnumber){
bool eventPass(std::vector<float>* jetPt, size_t jetnumber)
{
    float jetPtCut = 100.0;
    if(jetPt->size() < 1) return false;
    if(jetnumber > 200) return false;
    if((*jetPt)[0] < jetPtCut) return false;
    return true;
}

bool jetPass(std::vector< float > * jetEta, std::vector< float > * jetPt, int ijet){
    float jetPtCut = 100.0;
    float jetEtaCut = 1.6;
    if(fabs( (*jetEta)[ijet]) > jetEtaCut) return false;
    if((*jetPt)[ijet] < jetPtCut) return false;
    return true;
}