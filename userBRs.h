#ifndef USERBRS_H
#define USERBRS_H



#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "Quark.h"


Double_t tau0_to_U2(Double_t mN, Double_t tau0mN);

Double_t prodBR_lept(int idB, int idl, Double_t mN, Double_t tau0mN);
Double_t prodBR_semilept(int idB, int idl, int idH, Double_t mN, Double_t tau0mN);
Double_t decayBR_2body_semilept(int idl, int idH, Double_t mN, Double_t tau0mN);

#endif
