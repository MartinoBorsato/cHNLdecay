#ifndef USERBRS_H
#define USERBRS_H



#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "Quark.h"


//testing functions (error throwers)
void test_value(Double_t var, Double_t minval, Double_t maxval, std::string varname);

// utils
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::vector<double> > file_to_vec(std::string filename);
int getClosestIdx(Double_t val1, int idx1, Double_t val2, int idx2, Double_t target);
int findClosestIdx(std::vector<Double_t> arr, Double_t target);
Double_t tau0_to_U2(Double_t mN, Double_t tau0mN);
Double_t get_tau0ns(Double_t mN, Double_t U2);
Double_t get_GammaOverGammaInv(Double_t mN, Double_t U2);
Double_t get_GammaMesonOverGammaQuarks(Double_t mN, Double_t U2);

//BRs
Double_t prodBR_lept(int idB, int idl, Double_t mN, Double_t tau0mN);
Double_t prodBR_lept_U2(int idB, int idl, Double_t mN, Double_t U2);
Double_t prodBR_semilept(int idB, int idl, int idH, Double_t mN, Double_t tau0mN);
Double_t prodBR_semilept_U2(int idB, int idl, int idH, Double_t mN, Double_t U2);
//Double_t decayBR_2body_semilept(int idl, int idH, Double_t mN, Double_t tau0mN);
Double_t decayBR_lepton_meson(int idl, int idH, Double_t mN, Double_t tau0mN);
Double_t decayBR_lepton_lepton_neutrino(int idlA, int idlB, int idnuB, Double_t mN, Double_t tau0mN);

Double_t prodBR_Lambdab_proton(Double_t HNLmass, Double_t tau0);
Double_t prodBR_Lambdab_Lambdac(Double_t HNLmass, Double_t tau0);
#endif
