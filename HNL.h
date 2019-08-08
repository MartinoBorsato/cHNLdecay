//  cHNLdecay  --  calculate decay widths of Heavy Neutral Leptons
//  Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@posteo.de>
//
//  This file is part of cHNLdecay.
//
//  cHNLdecay is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  cHNLdecay is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef HNL_H
#define HNL_H

#include "Config.h"
#include "Lepton.h"
#include "Meson.h"
#include "TString.h"
#include <iostream>
#include <map>
#include <utility>

class HNL {
public:
  HNL(TString n, Double_t m, Double_t U2, const Lepton &a) {
    name = n;
    mass = m;
    angle = U2;
    generation.emplace_back(Lepton(a));
    //dc_c;
    majorana = true;
  }
  HNL(TString n, Double_t m, Double_t U2, const std::vector<Lepton> &a) {
    name = n;
    mass = m;
    angle = U2;
    for (auto l : a) {
      generation.emplace_back(Lepton(l));
    }
    //dc_c;

    majorana = true;
  }

  bool mixesWith(const Lepton &a) const {
    bool mixes = false;
    for (auto g : generation) {
      if (g == a)
        mixes = true;
    }
    return mixes;
  }

  Double_t getMass() const { return mass; }

  Double_t getAngle() const { return angle; }

  void setMass(Double_t m) {
    clearDecayChannels();
    mass = m;
  }

  void setAngle(Double_t a) {
    clearDecayChannels();
    angle = a;
  }

  TString getName() const { return name; }

  void newDecayChannel(std::vector<Int_t> particles, Double_t value) {
    dc_c.insert(std::pair<std::vector<Int_t>, Double_t>(particles, value));
  }

  bool existDecayChannel(std::vector<Int_t> p) const {
    return (dc_c.count(p) > 0);
  }

  void clearDecayChannels() { dc_c.clear(); }

  bool isMajorana() const { return majorana; }

  void setMajorana(bool val) { majorana = val; }

  std::map<std::vector<Int_t>, Double_t> getDecayChannels() const {
    return dc_c;
  }

  std::vector<Lepton> getGeneration() const { return generation; }

  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha,
                                const Quark &beta);
  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha,
                           const Quark &beta, const Quark &gamma);
  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha,
                           const Lepton &beta, bool invisible = true);
  Double_t getPartialWidthInv(std::shared_ptr<Config> cfg, const Lepton &alpha,
                              const Lepton &beta);
  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha,
                           const Meson &m);
  Double_t getTotalWidth(std::shared_ptr<Config> cfg,
                         const std::vector<Lepton> &leptons,
                         const std::vector<Meson> &particles);
  Double_t getTotalWidth(std::shared_ptr<Config> cfg,
                         const std::vector<Lepton> &leptons,
                         const std::vector<Quark> &particles);

private:
  TString name;
  Double_t mass;
  Double_t angle;
  std::vector<Lepton> generation;
  std::map<std::vector<Int_t>, Double_t> dc_c; // decay channels stored
  bool majorana;
};
#endif
