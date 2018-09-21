// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef PLOTS_H
#define PLOTS_H
#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "TString.h"
#include <iostream>

void plot_meson_pw(std::shared_ptr<Config> cfg, Lepton alpha,
                   std::vector<Meson> mesons, HNL N, TString output,
                   Int_t lowMass, Int_t highMass, Int_t stepsize = 10);
void plot_br(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
             std::vector<Meson> mesons, HNL N, TString output, Int_t lowMass,
             Int_t highMass, Int_t stepsize = 10);
void plot_qcd_coupling(TString output);
void plot_qcd_correction(TString output);
void plot_br_low(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
                 std::vector<Meson> mesons, HNL N, TString output,
                 Int_t lowMass, Int_t highMass, Int_t stepsize = 10);
#endif
