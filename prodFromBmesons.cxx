#include "formFactors.h"
#include "prodFromBmesons.h"
#include "Config.h"
#include "HNL.h"
#include "Lepton.h"
#include "Logger.h"
#include "Meson.h"
#include "auxfunctions.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

#include <chrono>
#include <thread>


// +++++++++++++++++ Tool functions +++++++++++++++++++++++++
Double_t kallen(Double_t a, Double_t b, Double_t c) {
	//std::cout<<"kallen " << pow(a,2) + pow(b,2) + pow(c,2) - 2*(a*b + b*c + c*a) << std::endl;
	
	return pow(a,2) + pow(b,2) + pow(c,2) - 2*(a*b + b*c + c*a);
}

Double_t Lambda(Double_t xi, Double_t yhp, Double_t yl, Double_t yN) {
	//std::cout<<"Lambda " << sqrt(kallen(1, pow(yhp, 2), xi)) * sqrt(kallen(xi, pow(yN, 2), pow(yl, 2))) << std::endl;
	//if (kallen(1, pow(yhp, 2), xi) < 0) std::cout << " _____kallen(1, pow(yhp, 2), xi) < 0_____" << kallen(1, pow(yhp, 2), xi) << std::endl;
	//if (kallen(xi, pow(yN, 2), pow(yl, 2)) < 0){
		 //std::cout << " _____kallen(xi, pow(yN, 2), pow(yl, 2)) < 0_____ = " << kallen(xi, pow(yN, 2), pow(yl, 2)) << std::endl;
		 //std::cout << " xi " << xi << std::endl;
		 //std::cout << " yhp " << yhp << std::endl;
		 //std::cout << " yl " << yl << std::endl;
		 //std::cout << " yN " << yN << std::endl;
	 //}
	return sqrt(abs(kallen(1, pow(yhp, 2), xi))) * sqrt(abs(kallen(xi, pow(yN, 2), pow(yl, 2))));
	//return kallen(1, pow(yhp, 2), xi) * kallen(xi, pow(yN, 2), pow(yl, 2));
}

Double_t Gminus(Double_t xi, Double_t yhp, Double_t yl, Double_t yN) {
	//std::cout<<"Lambda " << sqrt(kallen(1, pow(yhp, 2), xi)) * sqrt(kallen(xi, pow(yN, 2), pow(yl, 2))) << std::endl;
	return xi*(pow(yN, 2) +  pow(yl, 2)) - pow((pow(yN, 2) -  pow(yl, 2)),2) ;
}

Double_t Gplus(Double_t xi, Double_t yhp, Double_t yl, Double_t yN) {
	//std::cout<<"Lambda " << sqrt(kallen(1, pow(yhp, 2), xi)) * sqrt(kallen(xi, pow(yN, 2), pow(yl, 2))) << std::endl;
	return xi*(pow(yN, 2) +  pow(yl, 2)) + pow((pow(yN, 2) -  pow(yl, 2)),2) ;
}

Double_t F(Double_t xi, Double_t yhp) {
	//std::cout<<"Lambda " << sqrt(kallen(1, pow(yhp, 2), xi)) * sqrt(kallen(xi, pow(yN, 2), pow(yl, 2))) << std::endl;
	return pow(1-xi, 2) - 2*pow(yhp,2)*(1+xi) + pow(yhp,4);
}







// ++++++++++++++++ COMPUTATIONS FUNCTIONS +++++++++++++++++++++++

Double_t compute_fP1(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t xi_) {
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp = mesonp.getMass();
	Quark_Type u = meson.getU(); Quark_Type d = meson.getD();
	
	
	// Get quarks flavours to apply corresponding V_CKM
	
	Double_t mN = N.getMass();
	Double_t U2 = N.getAngle();
	
	Double_t q2_ = xi_ * pow(mh_,2);
	mpfr_t yhp, yl, yN, mh, q2, xi;
	// initialisation
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	mpfr_init2(q2, BITS);
	mpfr_init2(xi, BITS);
	
	// set values
	mpfr_set_d(yhp, mhp, MPFR_RNDD);
	mpfr_set_d(yl, ml, MPFR_RNDD);
	mpfr_set_d(yN, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_set_d(q2, q2_, MPFR_RNDD);
	mpfr_set_d(xi, xi_, MPFR_RNDD);
	
	mpfr_div(yhp, yhp, mh, MPFR_RNDD);
	mpfr_div(yl, yl, mh, MPFR_RNDD);
	mpfr_div(yN, yN, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD); // xi = q2 / mh2
	
	// Compute form factor for input q2
	mpfr_t fplus2;
	mpfr_init2(fplus2, BITS); 
	mpfr_set_d(fplus2, compute_ffactor(meson, mesonp, q2_, 1), MPFR_RNDD);
	mpfr_pow_ui(fplus2, fplus2, 2, MPFR_RNDD);
	
	//Compute Lambda3 term
	mpfr_t Lambda3;
	mpfr_init2(Lambda3, BITS); 
	Double_t Lambda_ = Lambda(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD));
	mpfr_set_d(Lambda3, Lambda_, MPFR_RNDD);
	mpfr_pow_ui(Lambda3, Lambda3, 3, MPFR_RNDD);
	//std::cout<<"__fPq : Gama3_" << pow(Lambda_,3) <<std::endl; 
	
	//Compute the whoe function, then divide it by xi 3
	mpfr_t xi3; 
	mpfr_init2(xi3, BITS); 
	mpfr_pow_ui(xi3, xi, 3, MPFR_RNDD);
	
	//std::cout << "xi: " << q2_/pow(mh_,2) << std::endl;
	//std::cout << "_xi3: " << pow(q2_/pow(mh_,2),3) << std::endl;
	//std::cout << "xi3: " << mpfr_get_d(xi3,MPFR_RNDD) << std::endl;
	
	//std::cout << "xi3: " << mpfr_get_d(xi3,MPFR_RNDD) << std::endl;
	
	mpfr_t res; 
	mpfr_init2(res, BITS); 
	
	//std::cout << "Lambda3: " << mpfr_get_d(Lambda3,MPFR_RNDD) << std::endl;
	//std::cout << "fplus2: " << mpfr_get_d(fplus2,MPFR_RNDD) << std::endl;
	//std::cout << "3*xi3: " << mpfr_get_d(xi3,MPFR_RNDD) << std::endl;
	
	mpfr_mul(res, fplus2, Lambda3, MPFR_RNDD);
	
	//std::cout << "res 1: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	mpfr_div(res, res, xi3, MPFR_RNDD); //////////////////////////////////
	mpfr_div_ui(res, res, 3, MPFR_RNDD);
	
	//std::cout << "fP1: res 2: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	//mpfr_div(res, res, 3., MPFR_RNDD);
	
	return mpfr_get_d(res, MPFR_RNDD);
									
}

Double_t compute_fP2(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t xi_) {
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp = mesonp.getMass();
	Quark_Type U = meson.getU(); Quark_Type D = meson.getD();
	
	// Get quarks flavours to apply corresponding V_CKM
	
	Double_t mN = N.getMass();
	Double_t U2 = N.getAngle();
	
	Double_t q2_ = xi_ * pow(mh_,2);
	
	mpfr_t yhp, yl, yN, mh, q2, xi;
	// initialisation
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	mpfr_init2(q2, BITS);
	mpfr_init2(xi, BITS);
	
	// set values
	mpfr_set_d(yhp, mhp, MPFR_RNDD);
	mpfr_set_d(yl, ml, MPFR_RNDD);
	mpfr_set_d(yN, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_set_d(q2, q2_, MPFR_RNDD);
	mpfr_set_d(xi, xi_, MPFR_RNDD);
	
	mpfr_div(yhp, yhp, mh, MPFR_RNDD);
	mpfr_div(yl, yl, mh, MPFR_RNDD);
	mpfr_div(yN, yN, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD); // xi = q2 / mh2
	
	// Compute form factor for input q2
	mpfr_t fplus2;
	mpfr_init2(fplus2, BITS); 
	mpfr_set_d(fplus2, compute_ffactor(meson, mesonp, q2_, 1), MPFR_RNDD);
	mpfr_pow_ui(fplus2, fplus2, 2, MPFR_RNDD);
	
	//Compute Lambda term
	mpfr_t Lambda_m;
	mpfr_init2(Lambda_m, BITS); 
	Double_t Lambda_ = Lambda(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD));
	mpfr_set_d(Lambda_m, Lambda_, MPFR_RNDD);
		//std::cout<<"__fP2 : Gama_" << Lambda_ <<std::endl; 
	//Compute Gminus term
	mpfr_t Gminus_m;
	mpfr_init2(Gminus_m, BITS); 
	Double_t Gminus_ = Gminus(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD));
	mpfr_set_d(Gminus_m, Gminus_, MPFR_RNDD);
		//std::cout<<"__fP2 : Gminus_" << Gminus_ <<std::endl;
	
	//Compute kallen(1,yhp2,xi) term
	mpfr_t kallen_m;
	mpfr_init2(kallen_m, BITS); 
	Double_t kallen_ = kallen(1., pow(mpfr_get_d(yhp, MPFR_RNDD),2), mpfr_get_d(xi, MPFR_RNDD));
	mpfr_set_d(kallen_m, kallen_, MPFR_RNDD);
		//std::cout<<"__fP2 : kallen_" << kallen_ <<std::endl; 
	
	//Compute the whole function, then divide it by 2*xi3
	mpfr_t xi3; 
	mpfr_init2(xi3, BITS); 
	mpfr_pow_ui(xi3, xi, 3, MPFR_RNDD);
	
	
	
	//std::cout << "xi3: " << mpfr_get_d(xi3,MPFR_RNDD) << std::endl;
	
	mpfr_t res; 
	mpfr_init2(res, BITS); 
	
	
	mpfr_mul(res, fplus2, Lambda_m, MPFR_RNDD);
		//std::cout << "res 1: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	mpfr_mul(res, res, Gminus_m, MPFR_RNDD);
		//std::cout << "res 2: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	mpfr_mul(res, res, kallen_m, MPFR_RNDD);
	
		//std::cout << "res 3: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	mpfr_div(res, res, xi3, MPFR_RNDD); /////////////////////////////
	mpfr_div_ui(res, res, 2, MPFR_RNDD);
	
	//std::cout << "fP2: res 4: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	//mpfr_div(res, res, 3., MPFR_RNDD);
	
	return mpfr_get_d(res, MPFR_RNDD);
									
}

Double_t compute_fP3(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t xi_) {
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp = mesonp.getMass();
	Quark_Type u = meson.getU(); Quark_Type d = meson.getD();
	
	// Get quarks flavours to apply corresponding V_CKM
	
	Double_t mN = N.getMass();
	Double_t U2 = N.getAngle();
	
	Double_t q2_ = xi_ * pow(mh_,2);
	
	mpfr_t yhp, yl, yN, mh, q2, xi;
	// initialisation
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	mpfr_init2(q2, BITS);
	mpfr_init2(xi, BITS);
	
	// set values
	mpfr_set_d(yhp, mhp, MPFR_RNDD);
	mpfr_set_d(yl, ml, MPFR_RNDD);
	mpfr_set_d(yN, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_set_d(q2, q2_, MPFR_RNDD);
	mpfr_set_d(xi, xi_, MPFR_RNDD);
	
	mpfr_div(yhp, yhp, mh, MPFR_RNDD);
	mpfr_div(yl, yl, mh, MPFR_RNDD);
	mpfr_div(yN, yN, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD);
	//mpfr_div(xi, xi, mh, MPFR_RNDD); // xi = q2 / mh2
	
	// Compute form factor for input q2
	mpfr_t fzero2;
	mpfr_init2(fzero2, BITS); 
	mpfr_set_d(fzero2, compute_ffactor(meson, mesonp, q2_, 0), MPFR_RNDD);
	mpfr_pow_ui(fzero2, fzero2, 2, MPFR_RNDD);
	
	//Compute Lambda term
	mpfr_t Lambda_m;
	mpfr_init2(Lambda_m, BITS); 
	Double_t Lambda_ = Lambda(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD));
	mpfr_set_d(Lambda_m, Lambda_, MPFR_RNDD);
	
	//Compute Gminus term
	mpfr_t Gminus_m;
	mpfr_init2(Gminus_m, BITS); 
	Double_t Gminus_ = Gminus(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD));
	mpfr_set_d(Gminus_m, Gminus_, MPFR_RNDD);
		//std::cout<<"__fP3 : Gminus_" << Gminus_ <<std::endl;
	
	//Compute (1-yhp^2)^2 term
	mpfr_t temp;
	mpfr_init2(temp, BITS); 
	mpfr_pow_ui(temp, yhp, 2, MPFR_RNDD);
	mpfr_ui_sub(temp, 1, temp, MPFR_RNDD);
	mpfr_pow_ui(temp, temp, 2, MPFR_RNDD);
	
	//Compute the whole function, then divide it by 2*xi3
	mpfr_t xi3; 
	mpfr_init2(xi3, BITS); 
	mpfr_pow_ui(xi3, xi, 3, MPFR_RNDD);
	
	
	
		//std::cout << "xi3: " << mpfr_get_d(xi3,MPFR_RNDD) << std::endl;
	
	mpfr_t res; 
	mpfr_init2(res, BITS); 
	
	
	mpfr_mul(res, fzero2, Lambda_m, MPFR_RNDD);
		//std::cout << "res 1: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	mpfr_mul(res, res, Gminus_m, MPFR_RNDD);
		//std::cout << "res 2: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	mpfr_mul(res, res, temp, MPFR_RNDD);
	
		//std::cout << "res 3: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	mpfr_div(res, res, xi3, MPFR_RNDD); ///////////////
	mpfr_div_ui(res, res, 2, MPFR_RNDD);
	
	//std::cout << "fP3: res 4: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	//mpfr_div(res, res, 3., MPFR_RNDD);
	
	return mpfr_get_d(res, MPFR_RNDD);
									
}


// +++++++++++++++++ INTEGRATION FUNCTION ++++++++++++++++++++++++++++++

Double_t integral_fP(Double_t a, Double_t b, Double_t nsteps, std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, int fct){
	
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t step = (b-a)/nsteps;
	if(step<0.){
		std::cerr<<"-------\nForbidden kinematics for mN = " << N.getMass() << " for the decay " << meson.getPdgId() << " to " << mesonp.getPdgId() <<" " << l.getPdgId()<<" N. Please change input HNL features.\n--------"<<std::endl;
		//	std::cout << "step: " << step << std::endl;
		//	std::cout << "mN: " << N.getMass() << std::endl;
			return 0.;
		}
	mpfr_t res;
	mpfr_init2(res, BITS); 
	mpfr_set_d(res, 0., MPFR_RNDD);
	
	//std::cout << "step: " << step << std::endl;
	
	//loop on nsteps steps and add integral using trapezoidal computation @each step
	Double_t x1, x2, x3; 
	Double_t fstep_, f1, f2, f3; 
	Double_t xi1, xi2, xi3;
	Double_t mh_ = meson.getMass();
	for(int i(0); i<nsteps; ++i){
		xi1 = a+i*step;
		xi2 = a+(i+0.5)*step;
		xi3 = a+(i+1.)*step;
		
		if(fct==1){
			f1=compute_fP1(cfg, N, l, meson, mesonp, xi1);
			f2=compute_fP1(cfg, N, l, meson, mesonp, xi2);
			f3=compute_fP1(cfg, N, l, meson, mesonp, xi3);
		}
		
		else if(fct==2){
			f1=compute_fP2(cfg, N, l, meson, mesonp, xi1);
			f2=compute_fP2(cfg, N, l, meson, mesonp, xi2);
			f3=compute_fP2(cfg, N, l, meson, mesonp, xi3);
		}
		
		else if(fct==3){
			f1=compute_fP3(cfg, N, l, meson, mesonp, xi1);
			f2=compute_fP3(cfg, N, l, meson, mesonp, xi2);
			f3=compute_fP3(cfg, N, l, meson, mesonp, xi3);
		}
		
		
		fstep_ = step*(1./6*f1 + 4./6*f2 + 1./6*f3);
		
		//std::cout << "ffstep_: " << fstep_ << std::endl;
		mpfr_t fstep;
		mpfr_init2(fstep, BITS); 
		mpfr_set_d(fstep, fstep_, MPFR_RNDD);
		
		mpfr_add(res, res, fstep, MPFR_RNDD);
	}
	
	return mpfr_get_d(res, MPFR_RNDD);
}



// +++++++++++++ LEPTONIC CASE +++++++++++++++++++++++++++

Double_t pw_prodFromBmeson_leptonic(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson){
	
	mpfr_t fermiC, fermiCsq, pi, fh;
	unsigned int BITS = cfg->getBITS();
	Double_t fh_ = meson.getDecayConstant(); 
	//std::cout<<"fh: " << fh_ << std::endl;
	mpfr_init2(fermiC, BITS);
	mpfr_init2(fermiCsq, BITS);
	mpfr_init2(pi, BITS);
	mpfr_init2(fh, BITS);
	mpfr_set_d(fh, fh_, MPFR_RNDD);
	cfg->getFermiC(fermiC);
	cfg->getPi(pi);
	
	/** Get high precision values **/
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass();
	Double_t mN = N.getMass();
	Double_t U2_ = N.getAngle();
	
	// Compute corresponding V_CKM
	
	
	//Double_t Vud_ = get_VCKM(Din,Dout); 
	//Double_t Vud_ = cfg->getVUDsq(Din, Dout); 
	Double_t Vud_ = 1.;
	//std::cout<<"meson name: " << meson.getName() << std::endl;
	//std::cout<<"mesonp name: " << mesonp.getName() << std::endl;
	//std::cout << "VCKM elmt: " << Vud_ << std::endl;
	
	mpfr_t yl2, yN2, mh, q2, xi, U2, Vud2;
	// initialisation
	
	mpfr_init2(yl2, BITS); 
	mpfr_init2(yN2, BITS);
	mpfr_init2(mh, BITS);
	//mpfr_init2(q2, BITS);
	//mpfr_init2(xi, BITS);
	mpfr_init2(U2, BITS);
	mpfr_init2(Vud2, BITS);
	
	// set values
	
	mpfr_set_d(yl2, ml, MPFR_RNDD);
	mpfr_set_d(yN2, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_div(yl2, yl2, mh, MPFR_RNDD);
	mpfr_div(yN2, yN2, mh, MPFR_RNDD);
	mpfr_pow_ui(yl2, yl2, 2, MPFR_RNDD);
	mpfr_pow_ui(yN2, yN2, 2, MPFR_RNDD);
	
//	std::cout<<"_______ VARIABLES _________" << std::endl;
//	std::cout<<"yl2: " << mpfr_get_d(yl2, MPFR_RNDD) << std::endl;
//	std::cout<<"yN2: " << mpfr_get_d(yN2, MPFR_RNDD) << std::endl;
	
	
	//mpfr_set_d(q2, q2_, MPFR_RNDD);
	//mpfr_set_d(xi, q2_, MPFR_RNDD);
	
	mpfr_set_d(U2, U2_, MPFR_RNDD);
	mpfr_set_d(Vud2, Vud_, MPFR_RNDD);
	
	/** Compute factor **/
	mpfr_t factor, tmp;
	mpfr_init2(factor, BITS);
	mpfr_init2(tmp, BITS);

	
	mpfr_pow_ui(factor, fermiC, 2, MPFR_RNDD);
//	std::cout<<"factor 1: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_pow_ui(tmp, mh, 3, MPFR_RNDD);
//	std::cout<<"factor 2: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_mul(factor, factor, tmp, MPFR_RNDD);
//	std::cout<<"factor 3: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_pow_ui(tmp, fh, 2, MPFR_RNDD);
//	std::cout<<"factor 4: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_mul(factor, factor, tmp, MPFR_RNDD);
//	std::cout<<"factor 5: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_mul(factor, factor, U2, MPFR_RNDD);
//	std::cout<<"factor 6: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_mul(factor, factor, Vud2, MPFR_RNDD);
//	std::cout<<"factor 7: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_div(factor, factor, pi, MPFR_RNDD);
//	std::cout<<"factor 8: " << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
	mpfr_div_ui(factor, factor, 8, MPFR_RNDD);

	
	/** Compute function(y_particles) **/
	
	mpfr_t factor2, tmp2, tmp22;
	mpfr_init2(factor2, BITS);
	mpfr_init2(tmp2, BITS);
	mpfr_init2(tmp22, BITS);
	
	mpfr_sub(tmp2, yN2, yl2, MPFR_RNDD);
	mpfr_pow_ui(tmp2, tmp2, 2, MPFR_RNDD);
	mpfr_add(tmp22, yN2, yl2, MPFR_RNDD);
	
	mpfr_sub(factor2, tmp22, tmp2, MPFR_RNDD);
	
	Double_t kln = sqrt(kallen(1,mpfr_get_d(yN2, MPFR_RNDD),mpfr_get_d(yl2, MPFR_RNDD)));

	
	mpfr_mul_d(factor2, factor2, kln,MPFR_RNDD);
	/** Compute whole PW value **/
	mpfr_mul(factor, factor, factor2, MPFR_RNDD);

	if(mpfr_get_d(factor, MPFR_RNDD)<0.) return 0.;
	else return mpfr_get_d(factor, MPFR_RNDD);
	
}

// ++++++++++++++++ SEMILEPTONIC VECTOR CASE +++++++++++++++++++++++++++

Double_t compute_fV(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t xi_, int part) {
	
	/** 
	 * 1: Vg2	2: Vf2	 3: Va+2	4: Va-2
	 * 5: Vfa+	6: Vfa-	 7: Va+a-
	**/
	
	unsigned int BITS = cfg->getBITS();
	
	// Masses in GeV
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp = mesonp.getMass();
	Quark_Type U = meson.getU(); Quark_Type D = meson.getD();
	
	// Get quarks flavours to apply corresponding V_CKM
	
	Double_t mN = N.getMass();//*1e-3;
	Double_t U2 = N.getAngle();
	
	Double_t q2_ = xi_ * pow(mh_,2);
	
	mpfr_t yhp, yl, yN, mh, q2, xi, prefactor, tmp1, tmp2, res1, res2, res;
	
	// initialisation
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	mpfr_init2(q2, BITS);
	mpfr_init2(xi, BITS);
	
	mpfr_init2(prefactor, BITS);
	mpfr_init2(tmp1, BITS);
	mpfr_init2(tmp2, BITS);
	mpfr_init2(res, BITS);
	mpfr_init2(res1, BITS);
	mpfr_init2(res2, BITS);
	
	
	// set values
	mpfr_set_d(yhp, mhp, MPFR_RNDD);
	mpfr_set_d(yl, ml, MPFR_RNDD);
	mpfr_set_d(yN, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_set_d(q2, q2_, MPFR_RNDD);
	mpfr_set_d(xi, xi_, MPFR_RNDD);
	
	mpfr_div(yhp, yhp, mh, MPFR_RNDD);
	mpfr_div(yl, yl, mh, MPFR_RNDD);
	mpfr_div(yN, yN, mh, MPFR_RNDD);
	
	mpfr_t xi2, xi3; 
	mpfr_init2(xi2, BITS); mpfr_init2(xi3, BITS); 
	mpfr_pow_ui(xi2, xi, 2, MPFR_RNDD); mpfr_pow_ui(xi3, xi, 3, MPFR_RNDD);
	//std::cout<<"xi3"<<mpfr_get_d(xi3, MPFR_RNDD)<<std::endl;
	
	// Compute form factors for input q2
	mpfr_t g, g2, f, f2, aplus, aminus, aplus2, aminus2;
	mpfr_init2(g, BITS); mpfr_init2(g2, BITS);
	mpfr_init2(f, BITS); mpfr_init2(f2, BITS);
	mpfr_init2(aplus, BITS); mpfr_init2(aplus2, BITS);
	mpfr_init2(aminus, BITS); mpfr_init2(aminus2, BITS);
	
	
	mpfr_set_d(g, compute_ffactor(meson, mesonp, q2_, 2), MPFR_RNDD);
	mpfr_set_d(f, compute_ffactor(meson, mesonp, q2_, 3), MPFR_RNDD);
	mpfr_set_d(aplus, compute_ffactor(meson, mesonp, q2_, 4), MPFR_RNDD);
	mpfr_set_d(aminus, compute_ffactor(meson, mesonp, q2_, 5), MPFR_RNDD);
	
	mpfr_pow_ui(g2, g, 2, MPFR_RNDD);
	mpfr_pow_ui(f2, f, 2, MPFR_RNDD);
	mpfr_pow_ui(aplus2, aplus, 2, MPFR_RNDD);
	mpfr_pow_ui(aminus2, aminus, 2, MPFR_RNDD);
			
			
	//Compute Lambda term
	mpfr_t Lambda_;
	mpfr_init2(Lambda_, BITS); 
	mpfr_set_d(Lambda_, Lambda(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD)), MPFR_RNDD);
	
	//Compute Gminus term
	mpfr_t Gminus_;
	mpfr_init2(Gminus_, BITS); 
	mpfr_set_d(Gminus_, Gminus(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD)), MPFR_RNDD);
	
	//Compute Gplus term
	mpfr_t Gplus_;
	mpfr_init2(Gplus_, BITS); 
	mpfr_set_d(Gplus_, Gplus(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD), mpfr_get_d(yl, MPFR_RNDD), mpfr_get_d(yN, MPFR_RNDD)), MPFR_RNDD);
	
	//Compute F term
	mpfr_t F_;
	mpfr_init2(F_, BITS); 
	//std::cout<< "F:" << F(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD)) <<std::endl;
	mpfr_set_d(F_, F(mpfr_get_d(xi, MPFR_RNDD), mpfr_get_d(yhp, MPFR_RNDD)), MPFR_RNDD);
			
	//Compute kallen(1,yhp2,xi) term
	mpfr_t kallen_;
	mpfr_init2(kallen_, BITS); 
	mpfr_set_d(kallen_, kallen(1., pow(mpfr_get_d(yhp, MPFR_RNDD),2), mpfr_get_d(xi, MPFR_RNDD)), MPFR_RNDD);
	
	
	switch(part){
		case 1:{
			// Prefactor
			mpfr_pow_ui(prefactor, mh, 2, MPFR_RNDD);
			mpfr_pow_ui(tmp1, yhp, 2, MPFR_RNDD);
			mpfr_mul(prefactor, prefactor, tmp1, MPFR_RNDD);
			mpfr_div_d(prefactor, prefactor, 3, MPFR_RNDD);
			//std::cout << "case 1a: " << mpfr_get_d(prefactor, MPFR_RNDD) << std::endl;
			
			// Inside integral
			mpfr_mul_ui(res1, xi2, 2, MPFR_RNDD);
			mpfr_sub(res1, res1, Gplus_, MPFR_RNDD);
			//std::cout << "case 1 a: " << mpfr_get_d(res1, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res1, F_, MPFR_RNDD);
			//std::cout << "case 1 b: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			//std::cout << "case 1 c: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res, g2, MPFR_RNDD);
			//std::cout << "case 1 d: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_div(res, res, xi2, MPFR_RNDD);
			//std::cout << "case 1 e: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			
			///\std::cout << "case 1: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		
		case 2:{
			// Prefactor
			mpfr_pow_ui(tmp1, mh, 2, MPFR_RNDD);
			mpfr_ui_div(prefactor, 1, tmp1, MPFR_RNDD);
			mpfr_div_ui(prefactor, prefactor, 24, MPFR_RNDD);
			//std::cout << "case 2a: " << mpfr_get_d(prefactor, MPFR_RNDD) << std::endl;
			// Inside integral
			mpfr_mul_ui(res1, xi2, 2, MPFR_RNDD);
			mpfr_sub(res1, res1, Gplus_, MPFR_RNDD);
			
			mpfr_pow_ui(res, yhp, 2, MPFR_RNDD);
			mpfr_mul(res, res, res1, MPFR_RNDD);
			mpfr_mul_ui(res, res, 12, MPFR_RNDD);
			mpfr_mul(res, res, xi, MPFR_RNDD);
			
			mpfr_pow_ui(res1, Lambda_, 2, MPFR_RNDD);
			mpfr_sub(res, res, res1, MPFR_RNDD);
			
			mpfr_pow_ui(res1, yl, 2, MPFR_RNDD);
			mpfr_pow_ui(res2, yN, 2, MPFR_RNDD);
			mpfr_sub(res1, res1, res2, MPFR_RNDD);
			mpfr_pow_ui(res1, res1, 2, MPFR_RNDD);
			mpfr_sub(res1, xi2, res1, MPFR_RNDD);
			mpfr_mul(res1, res1, F_, MPFR_RNDD);
			mpfr_mul_ui(res1, res1, 3, MPFR_RNDD);
			
			mpfr_add(res, res, res1, MPFR_RNDD);
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, f2, MPFR_RNDD);
			
			mpfr_div(res, res, xi3, MPFR_RNDD);
			
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			//std::cout << "case 2: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		
		case 3:{
			// Prefactor
			mpfr_pow_ui(prefactor, mh, 2, MPFR_RNDD);
			mpfr_div_ui(prefactor, prefactor, 24, MPFR_RNDD);
			//std::cout << "case 3a: " << mpfr_get_d(prefactor, MPFR_RNDD) << std::endl;
			// Inside integral
			mpfr_pow_ui(res, yhp, 2, MPFR_RNDD);
			mpfr_ui_sub(res, 1, res, MPFR_RNDD);
			mpfr_pow_ui(res, res, 2, MPFR_RNDD);
			//std::cout << "case 3b: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res, Gminus_, MPFR_RNDD);
			//std::cout << "case 3c: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul_ui(res, res, 3, MPFR_RNDD);
			//std::cout << "case 3d: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			
			mpfr_mul_ui(res1, xi2, 2, MPFR_RNDD);
			mpfr_sub(res1, res1, Gplus_, MPFR_RNDD);
			mpfr_mul(res1, res1, F_, MPFR_RNDD);
			
			mpfr_add(res, res, res1, MPFR_RNDD);
			
			
			mpfr_mul(res, res, F_, MPFR_RNDD);
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, aplus2, MPFR_RNDD);
			
			mpfr_div(res, res, xi3, MPFR_RNDD);
			
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			
			//std::cout << "case 3: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		case 4:{
			// Prefactor
			mpfr_pow_ui(prefactor, mh, 2, MPFR_RNDD);
			mpfr_div_d(prefactor, prefactor, 8, MPFR_RNDD);
			
			mpfr_mul(res, Gminus_, F_, MPFR_RNDD);
			//std::cout << "case 4a: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, aminus2, MPFR_RNDD);
			//mpfr_mul(res, res, f2, MPFR_RNDD);
			mpfr_div(res, res, xi, MPFR_RNDD);
			
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			//std::cout << "case 4: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		case 5:{
			// Prefactor
			mpfr_set_d(prefactor, 1, MPFR_RNDD);
			mpfr_div_d(prefactor, prefactor, 12, MPFR_RNDD);
			
			mpfr_pow_ui(res1, yl, 2, MPFR_RNDD);
			mpfr_pow_ui(res2, yN, 2, MPFR_RNDD);
			mpfr_sub(res1, res1, res2, MPFR_RNDD);
			mpfr_pow_ui(res1, res1, 2, MPFR_RNDD);
			
			mpfr_sub(res1, xi2, res1, MPFR_RNDD);
			mpfr_mul(res, res1, F_, MPFR_RNDD);
			mpfr_mul_ui(res, res, 3, MPFR_RNDD);
			mpfr_pow_ui(res1, Lambda_, 2, MPFR_RNDD);
			mpfr_sub(res, res, res1, MPFR_RNDD);
			
			mpfr_pow_ui(res1, yhp, 2, MPFR_RNDD);
			mpfr_ui_sub(res1, 1, res1, MPFR_RNDD);
			mpfr_sub(res1, res1, xi, MPFR_RNDD);
			
			mpfr_mul(res, res1, res, MPFR_RNDD);
			
			mpfr_mul(res1, Gminus_, F_, MPFR_RNDD);
			mpfr_mul(res1, res1, xi, MPFR_RNDD);
			mpfr_mul_ui(res1, res1, 3, MPFR_RNDD);
			
			mpfr_add(res, res, res1, MPFR_RNDD);
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, aplus, MPFR_RNDD);
			mpfr_mul(res, res, f, MPFR_RNDD);
			
			mpfr_div(res, res, xi3, MPFR_RNDD);
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			
			//std::cout << "case 5: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		case 6:{
			// Prefactor
			mpfr_set_d(prefactor, 1, MPFR_RNDD);
			mpfr_div_d(prefactor, prefactor, 4, MPFR_RNDD);
			
			//Inside integral
			mpfr_mul(res, Gminus_, F_, MPFR_RNDD);
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, aminus, MPFR_RNDD);
			mpfr_mul(res, res, f, MPFR_RNDD);
			mpfr_div(res, res, xi2, MPFR_RNDD);
			
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			//std::cout << "case 6: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		case 7:{
			// Prefactor
			mpfr_pow_ui(prefactor, mh, 2, MPFR_RNDD);
			mpfr_div_d(prefactor, prefactor, 4, MPFR_RNDD);
			
			mpfr_pow_ui(res, yhp, 2, MPFR_RNDD);
			mpfr_ui_sub(res, 1, yhp, MPFR_RNDD);
			mpfr_mul(res, res, Gminus_, MPFR_RNDD);
			mpfr_mul(res, res, F_, MPFR_RNDD);
			mpfr_mul(res, res, Lambda_, MPFR_RNDD);
			mpfr_mul(res, res, aminus, MPFR_RNDD);
			mpfr_mul(res, res, aplus, MPFR_RNDD);
			mpfr_div(res, res, xi2, MPFR_RNDD);
			
			mpfr_mul(res, res, prefactor, MPFR_RNDD);
			
			//std::cout << "case 7: " << mpfr_get_d(res, MPFR_RNDD) << std::endl;
			
			return mpfr_get_d(res, MPFR_RNDD);
			break;
		}
		
		default:
		{
			 std::cerr<<"Error: fV 'part' argument must be between 1 and 7" << std::endl; break;
		}
	}
	
	return 0.;
}


Double_t compute_integral(Double_t a, Double_t b, Double_t nsteps, std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp){
	if(b<a) return 0;
	else{
		Double_t step = (b-a)/nsteps;
		unsigned int BITS = cfg->getBITS();
		//vector<Double_t> fV_vec;
		Double_t xi1, xi2, xi3;
		Double_t fstep_, fV_xi1(0), fV_xi2(0), fV_xi3(0);
		mpfr_t res; 
		mpfr_init2(res, BITS); mpfr_set_d(res, 0, MPFR_RNDD); 
		mpfr_t fstep;
		mpfr_init2(fstep, BITS); 
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			for(size_t j(1); j<=7; ++j){
				fV_xi1 = fV_xi1 + compute_fV(cfg, N, l, meson, mesonp, xi1, j);
				//std::cout<<"fV_xi1, step " << j << ": " << fV_xi1 << std::endl;
				
				fV_xi2 = fV_xi2 + compute_fV(cfg, N, l, meson, mesonp, xi2, j);
				fV_xi3 = fV_xi3 + compute_fV(cfg, N, l, meson, mesonp, xi3, j);	
			}
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			//std::cout << "fstep_" << fstep_ << std::endl;
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			
			mpfr_add(res, res, fstep, MPFR_RNDD);
			
			
		}
		//std::cout << "res" << mpfr_get_d(res, MPFR_RNDD) << std::endl;
		
		return mpfr_get_d(res, MPFR_RNDD);
		//return 0;
	}
}

void MONITORING(Double_t nsteps, std::shared_ptr<Config> cfg, Lepton l, Meson meson, Meson mesonp){
	std::vector<double> mNvec = {0.,0.244, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25,  3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00,6.25,6.5};

	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp_ = mesonp.getMass();
	
		

	TString s = std::to_string(meson.getPdgId())+"_to_"+std::to_string(mesonp.getPdgId())+"_muN";
	std::ofstream ofile; 	ofile.open("../dat_files/integral_parts_"+s+".dat");

	Double_t step;
	
	unsigned int BITS = cfg->getBITS();
	//vector<Double_t> fV_vec;
	Double_t xi1, xi2, xi3;
	Double_t fstep_, fV_xi1(0), fV_xi2(0), fV_xi3(0);
	mpfr_t res; 
	mpfr_init2(res, BITS); mpfr_set_d(res, 0, MPFR_RNDD); 
	mpfr_t fstep;
	mpfr_init2(fstep, BITS); 
	
	Double_t mN;
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		ofile << mN<<" ";
	}
	ofile << std::endl;
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		
		if(step < 0) ofile << 0.0 <<" ";
		
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 1);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 1);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 1);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 2);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 2);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 2);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 3);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 3);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 3);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 4);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 4);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 4);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 5);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 5);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 5);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 6);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 6);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 6);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	for (int k(0); k<mNvec.size(); k++){
		mN = mNvec[k]*1e3;
		std::vector<Lepton> mixes_with; mixes_with={l};
		HNL N = HNL("HNL", mN, 1e-6, mixes_with);
		
		Double_t a = pow(ml/mh_ + mN/mh_,2); Double_t b = pow(1-(mhp_/mh_),2);
		step = (b-a)/nsteps;
		if(step < 0) ofile << 0.0 <<" ";
		else{
		for(int i(0); i<nsteps; ++i){
			xi1 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			xi2 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			xi3 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			fV_xi1=0; fV_xi2=0; fV_xi3=0;
			fstep_=0;
			
			fV_xi1 = compute_fV(cfg, N, l, meson, mesonp, xi1, 7);
			fV_xi2 = compute_fV(cfg, N, l, meson, mesonp, xi2, 7);
			fV_xi3 = compute_fV(cfg, N, l, meson, mesonp, xi3, 7);
			
			fstep_ = step*(1./6.*fV_xi1 + 4./6.*fV_xi2 + 1./6.*fV_xi3);
			
			mpfr_set_d(fstep, fstep_, MPFR_RNDD);
			mpfr_add(res, res, fstep, MPFR_RNDD);
		}
		ofile << mpfr_get_d(res, MPFR_RNDD)<<" ";
	}
	}
	ofile << std::endl;
	
	return;
}




Double_t pw_prodFromBmeson_semileptonic(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp){
	
	//std::cout << "0: pw_prodFromBmeson_semileptonic" << std::endl;
	
	bool process; // 0:pseudoscalar, 1:vector
	
	mpfr_t fermiC, fermiCsq, pi;
	unsigned int BITS = cfg->getBITS();
	mpfr_init2(fermiC, BITS);
	mpfr_init2(fermiCsq, BITS);
	mpfr_init2(pi, BITS);
	cfg->getFermiC(fermiC);
	cfg->getPi(pi);
	
	//std::cout << "1: pw_prodFromBmeson_semileptonic" << std::endl;
	//std::cout << "fermiC: " << mpfr_get_d(fermiC, MPFR_RNDD) << std::endl;
	
	/** Get high precision values **/
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp_ = mesonp.getMass();
	Double_t mN = N.getMass();
	Double_t U2_ = N.getAngle();
	
	// Compute corresponding V_CKM
	const Quark_Type Din = meson.getD(); const Quark_Type Dout = mesonp.getD();
	//std::cout << "2: pw_prodFromBmeson_semileptonic" << std::endl;
	//Double_t Vud_ = get_VCKM(Din,Dout); 
	//Double_t Vud_ = cfg->getVUDsq(Din, Dout); 
	Double_t Vud_ = 1.;
	//std::cout<<"meson name: " << meson.getName() << std::endl;
	//std::cout<<"mesonp name: " << mesonp.getName() << std::endl;
	//std::cout << "VCKM elmt: " << Vud_ << std::endl;
	
	mpfr_t yhp, yl, yN, mh, mhp, q2, xi, U2, Vud2;
	// initialisation
	
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	mpfr_init2(mhp, BITS);
	mpfr_init2(U2, BITS);
	mpfr_init2(Vud2, BITS);
	
	// set values
	
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	mpfr_set_d(mhp, mhp_, MPFR_RNDD);
	
	mpfr_d_div(yl, ml, mh, MPFR_RNDD);
	mpfr_d_div(yhp, mhp_, mh, MPFR_RNDD);
	mpfr_d_div(yN, mN, mh, MPFR_RNDD);
	
	//std::cout << "fermiC: " << mpfr_get_d(fermiC, MPFR_RNDD) << std::endl;
	//std::cout << "mh: " << mpfr_get_d(mh, MPFR_RNDD) << std::endl;
	//std::cout << "mhp: " << mpfr_get_d(mhp, MPFR_RNDD) << std::endl;
	
	//mpfr_set_d(q2, q2_, MPFR_RNDD);
	//mpfr_set_d(xi, q2_, MPFR_RNDD);
		
	mpfr_set_d(U2, U2_, MPFR_RNDD);
	mpfr_set_d(Vud2, Vud_, MPFR_RNDD);
	mpfr_pow_ui(Vud2, Vud2, 2, MPFR_RNDD);
		
	/** Compute factor **/
	mpfr_t factor, tmp;
	mpfr_init2(factor, BITS);
	mpfr_init2(tmp, BITS);
	
	// PSEUDOSCALAR
	if(mesonp.getMesonType()==MesonType::pseudoscalar){
	//if(1){	
		//std::cout << "MesonType::pseudoscalar" << std::endl;
		
		mpfr_pow_ui(factor, fermiC, 2, MPFR_RNDD);
		mpfr_pow_ui(tmp, mh, 5, MPFR_RNDD);
		mpfr_mul(factor, factor, tmp, MPFR_RNDD);
		mpfr_mul(factor, factor, U2, MPFR_RNDD);
		mpfr_mul(factor, factor, Vud2, MPFR_RNDD);
		mpfr_pow_ui(tmp, pi, 3, MPFR_RNDD);
		mpfr_div(factor, factor, tmp, MPFR_RNDD);
		mpfr_div_ui(factor, factor, 64, MPFR_RNDD);
			
		/** Compute 3 parts of the integral **/
		
		Double_t IP_(0);
		Double_t bmin = pow(ml/mh_ + mN/mh_,2); Double_t bmax = pow(1-(mhp_/mh_),2);
		for(int i(1); i<=3; ++i){
			IP_ += integral_fP(bmin, bmax, 5000, cfg, N, l, meson, mesonp, i);
			//std::cout << "INTEGRAL STEP " << i << " : " << IP_ << std::endl;
		}
		
		/** Compute whole PW value **/
		mpfr_t IP;
		mpfr_init2(IP, BITS); 
		mpfr_set_d(IP, IP_, MPFR_RNDD);
		
		mpfr_mul(factor, factor, IP, MPFR_RNDD);
		
		mpfr_clear(Vud2);

		if(mpfr_get_d(factor, MPFR_RNDD)<0.) return 0.;
		if(mesonp.getPdgId()==111) return 0.5*mpfr_get_d(factor, MPFR_RNDD);
		else return mpfr_get_d(factor, MPFR_RNDD);
	}
	
	// VECTOR
	else if(mesonp.getMesonType()==MesonType::vector){
	
		//std::cout<<"vector" << std::endl;
		mpfr_pow_ui(factor, fermiC, 2, MPFR_RNDD);
		//std::cout<<"fermic2" << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
		mpfr_pow_ui(tmp, mh, 7, MPFR_RNDD);
		mpfr_mul(factor, factor, tmp, MPFR_RNDD);
		mpfr_mul(factor, factor, U2, MPFR_RNDD);
		mpfr_mul(factor, factor, Vud2, MPFR_RNDD);
		mpfr_pow_ui(tmp, pi, 3, MPFR_RNDD);
		mpfr_div(factor, factor, tmp, MPFR_RNDD);
		mpfr_div_ui(factor, factor, 64, MPFR_RNDD);
		mpfr_pow_ui(tmp, mhp, 2, MPFR_RNDD);
		mpfr_div(factor, factor, tmp, MPFR_RNDD);
		//std::cout << "factor" << mpfr_get_d(factor, MPFR_RNDD) << std::endl;
		/** Compute 3 parts of the integral **/
		
		Double_t IP_(0);
		Double_t bmin = pow(ml/mh_ + mN/mh_,2); Double_t bmax = pow(1-(mhp_/mh_),2);
		
		IP_ = compute_integral(bmin, bmax, 500, cfg, N, l, meson, mesonp);
		
		//std::cout << "IP_: " << IP_ << std::endl;
		//IP_=1.;
		
		/** Compute whole PW value **/
		mpfr_t IP;
		mpfr_init2(IP, BITS); 
		mpfr_set_d(IP, IP_, MPFR_RNDD);
		//std::cout << "IP" << mpfr_get_d(IP, MPFR_RNDD) << std::endl;
		
		mpfr_mul(factor, factor, IP, MPFR_RNDD);
		
		mpfr_clears(fermiC, pi, Vud2, IP, yl, yN, yhp, mh, mhp, (mpfr_ptr)0);

		if(mpfr_get_d(factor, MPFR_RNDD)<0.) return 0.;
		if(mesonp.getPdgId()==213) return sqrt(0.5)*mpfr_get_d(factor, MPFR_RNDD);
		else return mpfr_get_d(factor, MPFR_RNDD);//*1.e-12;
	}
	
	else std::cerr<<"ERROR: MesonType not recognized (pseudosclar or vector)";
	return 0.;
	
}


// +++++++++++++ DISPLAY +++++++++++++++++++++++++++

void display_fform(HNL N, Lepton l, Meson meson, Meson mesonp){
	
	Double_t mh = 	meson.getMass();
	Double_t mhp = 	mesonp.getMass();
	Double_t ml = 	l.getMass();
	Double_t mN = 	N.getMass();
	Double_t q2;
	Double_t a = pow(ml/mh + mN/mh,2); Double_t b = pow(1-(mhp/mh),2);
	
	TString s = std::to_string(meson.getPdgId())+"_to_"+std::to_string(mesonp.getPdgId())+"_muN";
	std::ofstream ofile; 	ofile.open("../dat_files/ffactors_"+s+".dat");

	
	//Sequence the intervall and compute f
	int nsteps = 10;
	Double_t step = (b-a)/nsteps;
	
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << q2<<" ";
	}
	ofile << std::endl;
	// fplus
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,1)<<" ";
	}
	ofile << std::endl;
	// fzero
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	// g
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	// f
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	// aplus
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	// aminus
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	
	
	return;
}

///[1705.10606] and [1503.01421]////////////////////////////////////////

// hyperon decay (2 poss)
Double_t compute_dint_Lambdab_proton(std::shared_ptr<Config> cfg, HNL N, Lepton l, Double_t mh, Double_t mhp, Double_t q2) {
	
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t ml = l.getMass(); 
	Double_t mN = N.getMass();//*1e-3;
	Double_t U2 = N.getAngle();
	
	
	Double_t factor1; 
	
	// Form factors;
	
	Double_t f1V, f1A, f2V, f2A, f3V, f3A;
	Double_t a1V, a1A, a2V, a2A, a3V, a3A, a12V, a12A, a13V, a13A;
	
	f1V = f1V_Lambdab_proton(q2);	f1A = f1A_Lambdab_proton(q2);
	f2V = f2V_Lambdab_proton(q2);	f2A = f2A_Lambdab_proton(q2);
	f3V = f3V_Lambdab_proton(q2);	f3A = f3A_Lambdab_proton(q2);
	
	
	a1V = a1(q2, mh, mhp, ml, mN, 1);    a1A = a1(q2, mh, mhp, ml, mN, 0);
	a2V = a2(q2, mh, mhp, ml, mN, 1);    a2A = a2(q2, mh, mhp, ml, mN, 0);
	a3V = a3(q2, mh, mhp, ml, mN, 1);    a3A = a3(q2, mh, mhp, ml, mN, 0);
	a12V = a12(q2, mh, mhp, ml, mN, 1); a12A = a12(q2, mh, mhp, ml, mN, 0);
	a13V = a13(q2, mh, mhp, ml, mN, 1); a13A = a13(q2, mh, mhp, ml, mN, 0);
	

	factor1 = sqrt(kallen(pow(ml,2), pow(mN,2), q2)*kallen(pow(mh,2), pow(mhp,2), q2));
	
	return factor1*(
			 16/(3*pow(q2,3))*pow(f1V,2)*a1V
			+8/(3*pow(mh,2)*pow(q2,2))*pow(f2V,2)*a2V
			+8/(3*pow(mh,2)*q2)*pow(f3V,2)*a3V
			+32/(mh*pow(q2,2))*(
				 f1V*f2V*a12V
				+f1V*f3V*a13V)
			+16/(3*pow(q2,3))*pow(f1A,2)*a1A
			+8/(3*pow(mh,2)*pow(q2,2))*pow(f2A,2)*a2A
			+8/(3*pow(mh,2)*q2)*pow(f3A,2)*a3A
			+32/(mh*pow(q2,2))*(
				 f1A*f2A*a12A 
				+f1A*f3A*a13A
				)
			);
}

Double_t compute_dint_Lambdab_Lambdac(std::shared_ptr<Config> cfg, HNL N, Lepton l, Double_t mh, Double_t mhp, Double_t q2) {
	
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t ml = l.getMass(); 
	Double_t mN = N.getMass();//*1e-3;
	Double_t U2 = N.getAngle();
	
	
	Double_t factor1; 
	
	// Form factors;
	
	Double_t f1V, f1A, f2V, f2A, f3V, f3A;
	Double_t a1V, a1A, a2V, a2A, a3V, a3A, a12V, a12A, a13V, a13A;
	
	f1V = f1V_Lambdab_Lambdac(q2);	f1A = f1A_Lambdab_Lambdac(q2);
	f2V = f2V_Lambdab_Lambdac(q2);	f2A = f2A_Lambdab_Lambdac(q2);
	f3V = f3V_Lambdab_Lambdac(q2);	f3A = f3A_Lambdab_Lambdac(q2);
	
	
	a1V = a1(q2, mh, mhp, ml, mN, 1);    a1A = a1(q2, mh, mhp, ml, mN, 0);
	a2V = a2(q2, mh, mhp, ml, mN, 1);    a2A = a2(q2, mh, mhp, ml, mN, 0);
	a3V = a3(q2, mh, mhp, ml, mN, 1);    a3A = a3(q2, mh, mhp, ml, mN, 0);
	a12V = a12(q2, mh, mhp, ml, mN, 1); a12A = a12(q2, mh, mhp, ml, mN, 0);
	a13V = a13(q2, mh, mhp, ml, mN, 1); a13A = a13(q2, mh, mhp, ml, mN, 0);
	

	factor1 = sqrt(kallen(pow(ml,2), pow(mN,2), q2)*kallen(pow(mh,2), pow(mhp,2), q2));
	
	return factor1*(
			 16/(3*pow(q2,3))*pow(f1V,2)*a1V
			+8/(3*pow(mh,2)*pow(q2,2))*pow(f2V,2)*a2V
			+8/(3*pow(mh,2)*q2)*pow(f3V,2)*a3V
			+32/(mh*pow(q2,2))*(
				 f1V*f2V*a12V
				+f1V*f3V*a13V)
			+16/(3*pow(q2,3))*pow(f1A,2)*a1A
			+8/(3*pow(mh,2)*pow(q2,2))*pow(f2A,2)*a2A
			+8/(3*pow(mh,2)*q2)*pow(f3A,2)*a3A
			+32/(mh*pow(q2,2))*(
				 f1A*f2A*a12A 
				+f1A*f3A*a13A
				)
			);
	
}

///////////////////////////////

Double_t pw_prod_Lambdab_proton(std::shared_ptr<Config> cfg, HNL N, Lepton l){
	
	Double_t GF = 1e-11;
	Double_t pi = 3.141592;
	Double_t hbar = 1.;
	Double_t VCKM(1.);
	Double_t tauLb(1e-12);
	
	Double_t ml = l.getMass();
	Double_t mN = N.getMass();
	Double_t U2 = N.getAngle();
	
	Double_t mh = 5.;
	Double_t mhp = 1.;
	
	// Integration
	Double_t a = pow(ml+mN,2);
	Double_t b = pow(mh-mhp,2);
	Double_t nsteps = 1000;
	Double_t step = (b-a)/nsteps;
	Double_t q21, q22, q23;
	Double_t f_q21, f_q22, f_q23;
	Double_t fstep(0), res(0);
	
	if(step < 0) return 0;
	else{
		
		for(int i(0); i<nsteps; ++i){
			q21 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			q22 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			q23 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			
			f_q21 = compute_dint_Lambdab_proton(cfg, N, l, mh, mhp, q21);
			f_q22 = compute_dint_Lambdab_proton(cfg, N, l, mh, mhp, q22);
			f_q23 = compute_dint_Lambdab_proton(cfg, N, l, mh, mhp, q23);
			
			
			fstep = step*(1./6.*f_q21 + 4./6.*f_q22 + 1./6.*f_q23);
			
			res += fstep;
		}
		return pow(GF,2)*tauLb/(512*pow(pi,3)*pow(mh,3)*hbar)*pow(VCKM,2)*U2*res;
	}
}

Double_t pw_prod_Lambdab_Lambdac(std::shared_ptr<Config> cfg, HNL N, Lepton l){
	
	Double_t GF = 1e-11;
	Double_t pi = 3.141592;
	Double_t hbar = 1.;
	Double_t VCKM(1.);
	Double_t tauLb(1e-12);
	
	Double_t ml = l.getMass();
	Double_t mN = N.getMass();
	Double_t U2 = N.getAngle();
	
	Double_t mh = 5.;
	Double_t mhp = 1.;
	
	// Integration
	Double_t a = pow(ml+mN,2);
	Double_t b = pow(mh-mhp,2);
	Double_t nsteps = 1000;
	Double_t step = (b-a)/nsteps;
	Double_t q21, q22, q23;
	Double_t f_q21, f_q22, f_q23;
	Double_t fstep(0), res(0);
	
	if(step < 0) return 0;
	else{
		
		for(int i(0); i<nsteps; ++i){
			q21 = a+i*step; //std::cout << "xi1" << xi1 << std::endl;
			q22 = a+(i+0.5)*step; //std::cout << "xi2" << xi2 << std::endl;
			q23 = a+(i+1.)*step; //std::cout << "xi3" << xi3 << std::endl;
			
			
			f_q21 = compute_dint_Lambdab_Lambdac(cfg, N, l, mh, mhp, q21);
			f_q22 = compute_dint_Lambdab_Lambdac(cfg, N, l, mh, mhp, q22);
			f_q23 = compute_dint_Lambdab_Lambdac(cfg, N, l, mh, mhp, q23);
			

			fstep = step*(1./6.*f_q21 + 4./6.*f_q22 + 1./6.*f_q23);
			
			res += fstep;
		}
		return pow(GF,2)*tauLb/(512*pow(pi,3)*pow(mh,3)*hbar)*pow(VCKM,2)*U2*res;
	}
}




