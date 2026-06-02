#pragma once

#include "config.hpp"
#include "Math/Boost.h"

float GetMass(int pid){
    switch (abs(pid)){
        case 211:   return 0.13957;
        case 321:   return 0.49367;
        case 2212:  return 0.93827;
        case 11:    return 0.000511;
        case 13:    return 0.10566;
        case 22:    return 0.0;
        default:    return 0.13957;
    }
}

float GetQ(const ROOT::Math::PtEtaPhiMVector &p1, const ROOT::Math::PtEtaPhiMVector &p2){
    ROOT::Math::PtEtaPhiMVector diff = p1 - p2;
    double q = -(diff.M2());
    return (q >= 0.0 ? (float)TMath::Sqrt(q) : 0.f);
}

float GetQ_JetFrame(const ROOT::Math::PtEtaPhiMVector &p1,
                    const ROOT::Math::PtEtaPhiMVector &p2,
                    const ROOT::Math::PtEtaPhiMVector &jet)
{
    auto beta = jet.BoostToCM();
    ROOT::Math::Boost booster(beta.X(), beta.Y(), beta.Z());
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> lv1(p1);
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> lv2(p2);
    auto lv1b = booster(lv1);
    auto lv2b = booster(lv2);
    auto lvdiff = lv1b - lv2b;
    double q = -(lvdiff.M2());
    return (q >= 0.0 ? (float)TMath::Sqrt(q) : 0.f);
}

float GetQ_JetFrame_Double(const ROOT::Math::PtEtaPhiMVector &p1,
                           const ROOT::Math::PtEtaPhiMVector &jet1,
                           const ROOT::Math::PtEtaPhiMVector &p2,
                           const ROOT::Math::PtEtaPhiMVector &jet2)
{
    auto beta1 = jet1.BoostToCM();
    ROOT::Math::Boost booster1(beta1.X(), beta1.Y(), beta1.Z());
    auto beta2 = jet2.BoostToCM();
    ROOT::Math::Boost booster2(beta2.X(), beta2.Y(), beta2.Z());
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> lv1(p1);
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> lv2(p2);
    auto lv1b = booster1(lv1);
    auto lv2b = booster2(lv2);
    auto lvdiff = lv1b - lv2b;
    double q = -(lvdiff.M2());
    return (q >= 0.0 ? (float)TMath::Sqrt(q) : 0.f);
}