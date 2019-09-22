#ifndef PRODFROMBMESONS_H
#define PRODFROMBMESONS_H

#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "auxfunctions.h"
#include <vector>

Double_t compute_ffactor(Meson m, Meson mp, Double_t q2, int opt);



// hyperons decays

std::vector<Double_t> acoeffs_Lambdab_proton();
std::vector<Double_t> acoeffs_Lambdab_Lambdac();

Double_t f_Lambdab_proton(Double_t q2, int idx);  // f[idx] = fplus, f0, fperp, gplus, g0, gperp
Double_t f_Lambdab_Lambdac(Double_t q2, int idx); // f[idx] = fplus, f0, fperp, gplus, g0, gperp


Double_t f1V_Lambdab_proton(Double_t q2);
Double_t f2V_Lambdab_proton(Double_t q2);
Double_t f3V_Lambdab_proton(Double_t q2); 
Double_t f1A_Lambdab_proton(Double_t q2);
Double_t f2A_Lambdab_proton(Double_t q2);
Double_t f3A_Lambdab_proton(Double_t q2); 

Double_t f1V_Lambdab_Lambdac(Double_t q2);
Double_t f2V_Lambdab_Lambdac(Double_t q2);
Double_t f3V_Lambdab_Lambdac(Double_t q2); 
Double_t f1A_Lambdab_Lambdac(Double_t q2);
Double_t f2A_Lambdab_Lambdac(Double_t q2);
Double_t f3A_Lambdab_Lambdac(Double_t q2); 

Double_t Sigma_hyp(int sgn, Double_t mh, Double_t mhp);
Double_t Delta_hyp(int sgn, Double_t mh, Double_t mhp);

Double_t a1(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V);
Double_t a2(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V);
Double_t a3(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V);
Double_t a12(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V);
Double_t a13(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V);


#endif
