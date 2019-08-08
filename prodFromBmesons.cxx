
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
	if (kallen(1, pow(yhp, 2), xi) < 0) std::cout << " _____kallen(1, pow(yhp, 2), xi) < 0_____" << kallen(1, pow(yhp, 2), xi) << std::endl;
	if (kallen(xi, pow(yN, 2), pow(yl, 2)) < 0){
		 std::cout << " _____kallen(xi, pow(yN, 2), pow(yl, 2)) < 0_____ = " << kallen(xi, pow(yN, 2), pow(yl, 2)) << std::endl;
		 std::cout << " xi " << xi << std::endl;
		 std::cout << " yhp " << yhp << std::endl;
		 std::cout << " yl " << yl << std::endl;
		 std::cout << " yN " << yN << std::endl;
	 }
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


// ------ Meson form factors
Double_t compute_ffactor(Meson h, Meson hp, Double_t q2, bool opt) {
	
	// Uses equation (C.16) of arXiv paper
	
	Double_t 	Mpole, mh, mhp, a0, a1, a2, z, tplus, t0, factor, sum;
	Int_t 		IDh, IDhp;
	
	mh =  h.getMass(); 		mhp =  hp.getMass();
	IDh = h.getPdgId();		IDhp = hp.getPdgId();
	
	if(opt){		//form factor fplus
				
		switch(IDhp){	
			case 411:			//D
			case 413:			//Dstar
			case 421:			//D0bar
			case 423:			//D0starbar
			case 431:			//Ds
			case 433:			//Dsstar
				//Mpole = inf;
				a0 = 0.909;
				a1 = -7.11;
				a2 = 66;
				factor = 1.;
			case 313:			//K
			case 323:			//Kstar
				Mpole = 5325;
				a0 = 0.360;
				a1 = -0.828;
				a2 = 1.1;
				factor = 1/(1-q2/pow(Mpole, 2));
			case 211:			//pi
			case 111:			//pi0
				Mpole = 5325;
				a0 = 0.404;
				a1 = -0.68;
				a2 = -0.86;	
				factor = 1/(1-q2/pow(Mpole, 2));
		}
	}
	
	else{						//form factor f0
		
		switch(IDhp){	
			case 411:			//D
			case 413:			//Dstar
			case 421:			//D0bar
			case 423:			//D0starbar
			case 431:			//Ds
			case 433:			//Dsstar
				//Mpole = inf;
				a0 = 0.794;
				a1 = -2.45;
				a2 = 33;
				factor = 1.;
			case 313:			//K
			case 323:			//Kstar
				Mpole = 5650;
				a0 = 0.233;
				a1 = 0.197;
				a2 = 0.18;
				factor = 1/(1-q2/pow(Mpole, 2));
				//std::cout<<"factor " <<std::endl;
			case 211:			//pi
			case 111:			//pi0
				Mpole = 5650;
				a0 = 0.490;
				a1 = -1.61;
				a2 = 0.93;	
				factor = 1/(1-q2/pow(Mpole, 2));
				
				//std::cout<<"factor " <<std::endl;
		}
	}
	
	tplus = pow(mh+mhp, 2);
	t0 	  = (mh+mhp)*pow(sqrt(mh)-sqrt(mhp),2);
	
	z = sqrt(tplus-q2) - sqrt(tplus-t0);
	z = z/(sqrt(tplus-q2) + sqrt(tplus-t0));
	
	sum  = a0;
	sum += a1*(z 		- 1/3*pow(z,3));
	sum += a2*(pow(z,2)	+ 2/3*pow(z,3));
	
	
	return factor*sum;
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
	
	//std::cout << "res 2: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
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
	
		//std::cout << "res 4: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
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
	
		//std::cout << "res 4: " << mpfr_get_d(res,MPFR_RNDD) << std::endl;
	
	//mpfr_div(res, res, 3., MPFR_RNDD);
	
	return mpfr_get_d(res, MPFR_RNDD);
									
}


// +++++++++++++++++ INTEGRATION FUNCTION ++++++++++++++++++++++++++++++

Double_t integral_fP(Double_t a, Double_t b, Double_t nsteps, std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, int fct){
	
	
	unsigned int BITS = cfg->getBITS();
	
	Double_t step = (b-a)/nsteps;
	if(step<0.){
			std::cout << "step: " << step << std::endl;
			std::cout << "mN: " << N.getMass() << std::endl;
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
		
		mpfr_t fstep;
		mpfr_init2(fstep, BITS); 
		mpfr_set_d(fstep, fstep_, MPFR_RNDD);
		
		mpfr_add(res, res, fstep, MPFR_RNDD);
	}
	
	return mpfr_get_d(res, MPFR_RNDD);
}


// +++++++++++++++++++++++++ PARTIAL WIDTH ++++++++++++++++++++++++++++++

Double_t pw_prodFromBmeson_semileptonic(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp){
	
	mpfr_t fermiC, fermiCsq, pi;
	unsigned int BITS = cfg->getBITS();
	mpfr_init2(fermiC, BITS);
	mpfr_init2(fermiCsq, BITS);
	mpfr_init2(pi, BITS);
	cfg->getFermiC(fermiC);
	cfg->getPi(pi);
	
	//std::cout << "fermiC: " << mpfr_get_d(fermiC, MPFR_RNDD) << std::endl;
	
	/** Get high precision values **/
	Double_t ml = l.getMass(); 
	Double_t mh_ = meson.getMass(); Double_t mhp = mesonp.getMass();
	Double_t mN = N.getMass();
	Double_t U2_ = N.getAngle();
	
	// Compute corresponding V_CKM
	const Quark_Type Din = meson.getD(); const Quark_Type Dout = mesonp.getD();
	
	//Double_t Vud_ = get_VCKM(Din,Dout); 
	//Double_t Vud_ = cfg->getVUDsq(Din, Dout); 
	Double_t Vud_ = 1.;
	//std::cout<<"meson name: " << meson.getName() << std::endl;
	//std::cout<<"mesonp name: " << mesonp.getName() << std::endl;
	//std::cout << "VCKM elmt: " << Vud_ << std::endl;
	
	mpfr_t yhp, yl, yN, mh, q2, xi, U2, Vud2;
	// initialisation
	mpfr_init2(yhp, BITS); 
	mpfr_init2(yl, BITS); 
	mpfr_init2(yN, BITS);
	mpfr_init2(mh, BITS);
	//mpfr_init2(q2, BITS);
	//mpfr_init2(xi, BITS);
	mpfr_init2(U2, BITS);
	mpfr_init2(Vud2, BITS);
	
	// set values
	mpfr_set_d(yhp, mhp, MPFR_RNDD);
	mpfr_set_d(yl, ml, MPFR_RNDD);
	mpfr_set_d(yN, mN, MPFR_RNDD);
	mpfr_set_d(mh, mh_, MPFR_RNDD);
	
	//mpfr_set_d(q2, q2_, MPFR_RNDD);
	//mpfr_set_d(xi, q2_, MPFR_RNDD);
	
	mpfr_set_d(U2, U2_, MPFR_RNDD);
	mpfr_set_d(Vud2, Vud_, MPFR_RNDD);
	mpfr_pow_ui(Vud2, Vud2, 2, MPFR_RNDD);
	
	/** Compute factor **/
	mpfr_t factor, tmp;
	mpfr_init2(factor, BITS);
	mpfr_init2(tmp, BITS);
	
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
	Double_t bmin = pow(ml/mh_ + mN/mh_,2); Double_t bmax = pow(1-(mhp/mh_),2);
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
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,1)<<" ";
	}
	ofile << std::endl;
	for(int i(0); i<nsteps; i++){
		q2 = pow(mh,2)*(a+step*i);
		ofile << compute_ffactor(meson,mesonp,q2,0)<<" ";
	}
	ofile << std::endl;
	
	
	
	return;
}

