#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "Quark.h"

#include "Config.h"
#include "Logger.h"
#include "ParticleCatalogue.h"
#include "auxfunctions.h"
#include "partialWidths.h"
#include "prodFromBmesons.h"
#include "userBRs.h"

#include "plots.h"
#include <fstream>
#include <getopt.h>
#include <gmp.h>
#include <iomanip>
#include <iostream>
#include <mpfr.h>
#include <stdlib.h>


// Global variables
auto cfg = std::make_shared<Config>(); // set BITS and initialize constants

double hbar=6.582119514e-10; //in MeV/s
std::vector<std::vector<Double_t>> VCKM{{0.97427, 0.22534, 0.00351}, 
										{0.22520, 0.97344, 0.04120}, 
										{0.00867, 0.04040, 0.999146}};

Double_t mN = 0;
Double_t U2 = 0;

bool majorana = true;

// variables used for searching for angle / parsing file for angles
bool isLoad = false;
TString loadPath;
Double_t ctau = 0;

// sensible default for loglevel
//Level gLOGLEVEL = Level::INFO;

ParticleCatalogue pc;

std::vector<Meson> mesons = pc.getAllMesons();
	

///////////////////////////////////////////////////////////////////
// Declaration of all the particles we need for the HNL production
	
	
// Leptons

Double_t muonMass = 105.6583715;      // MeV
Double_t electronMass = 0.5109989461; // MeV
Double_t tauMass = 1776.82;           // MeV

Lepton el = Lepton(11, electronMass);
Lepton mu = Lepton(13, muonMass);
Lepton tau = Lepton(15, tauMass);	

std::vector<Lepton> mixes_with={mu};
	
// Mothers mesons

//Double_t tau0B0, tau0Bp, tau0Bs, tau0Bc;

Double_t tau0B0 = 1.520e-12;		//s
Double_t tau0Bp = 1.638e-12;		//s
Double_t tau0Bs = 1.510e-12;		//s
Double_t tau0Bc = 0.507e-12;		//s


Meson Bp = Meson(521, 5279.29, 187, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::up, Quark_Type::bottom);
Meson B0 = Meson(511, 5279.61, 130, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::bottom, Quark_Type::down);
Meson Bs = Meson(531, 5366.79, 228, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::strange, Quark_Type::bottom);
Meson Bc = Meson(541, 6275.1, 434, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::charm, Quark_Type::bottom); 

// HNL


// Daughters mesons
	
Meson pi = 		Meson(211, 139.57018, 130.2, MesonType::pseudoscalar, Charge::charged, 
				Quark_Type::up, Quark_Type::down);
Meson pi0 = 	Meson(111, 139.57018, 130.2, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::up, Quark_Type::up);
Meson rho = 	Meson(213, 775.4, 208.5, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::up, Quark_Type::down);
Meson rho0 = 	Meson(113, 775.49, 208.5, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::up, Quark_Type::up);
Meson K = 		Meson(313, 493.677, 159.8, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::strange, Quark_Type::up);
Meson Kst = 	Meson(323, 891.92, 212, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::charm, Quark_Type::up);
Meson D = 		Meson(411, 1869.62, 222.6, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::down, Quark_Type::charm);
Meson Dstar = 	Meson(413, 2010.26, 535000, MesonType::pseudoscalar, Charge::charged, 
				Quark_Type::charm, Quark_Type::down);
Meson D0bar = 	Meson(421, 1864.84, 212, MesonType::pseudoscalar, Charge::charged,
				Quark_Type::up, Quark_Type::charm);
Meson D0barstar = Meson(423, 2006.97, 212, MesonType::pseudoscalar, Charge::neutral,
				Quark_Type::down, Quark_Type::charm);
Meson Ds = 		Meson(431, 1968.47, 249, MesonType::pseudoscalar, Charge::charged, 
				Quark_Type::charm, Quark_Type::strange);              
Meson Dstars = 	Meson(433, 2112.1, 650000, MesonType::vector, Charge::charged, 
				Quark_Type::charm, Quark_Type::strange);                 

//HNL N("tmp", 0.0, 0.0, mixes_with);        

/** FUNCTIONS **/

Double_t tau0_to_U2(Double_t mN, Double_t tau0mN){
	Double_t U2;
	//computations
	U2=1.;
	return U2;
}


// production BR: leptonic case

Double_t prodBR_lept(int idB, int idl, Double_t mN, Double_t tau0mN){
	
	Double_t BR;
	Double_t pw, totw, tau0B;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	//Declare the HNL
	HNL N = HNL("HNL", mN, U2, mixes_with);
	N.setMajorana(majorana);
	
	Meson B; Lepton l;
	int i,j;
	
	switch(idB){
		case 521: B = Bp; i = 2; j = 0; tau0B = tau0Bp; break;
		case 541: B = Bc; i = 2; j = 1; tau0B = tau0Bc; break;	
	}
		
	switch(idB){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
			
	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_leptonic(cfg, N, l, B);
	totw = hbar/tau0B;
	return pw/totw;
}

// production BR: semileptonic case in pseudoscalar meson 

Double_t prodBR_semilept(int idB, int idl, int idH, Double_t mN, Double_t tau0mN){
	
	Double_t BR;
	Double_t pw, totw, tau0B;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	//Declare the HNL
	HNL N = HNL("HNL", mN, U2, mixes_with);
	bool majorana = true;
	N.setMajorana(majorana);
	
	Meson B; Meson H; Lepton l;
	bool isP(1);
	int i,j;
	
	switch(idB){
		case 511: B = B0; j = 0; tau0B = tau0B0; break;
		case 521: B = Bp; i = 0; tau0B = tau0Bp; break;
		case 531: B = Bs; j = 1; tau0B = tau0Bs; break;
		case 541: B = Bc; i = 1; tau0B = tau0Bc; break;	
	}
	
	// try to do it with jkey
	switch(idH){
		case 111: H = pi; j = 0; break;

	}
		
	switch(idB){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
	
	if(isP)	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_semileptonic(cfg, N, l, B, H); // pseudoscalar meson
	else 	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_semileptonic(cfg, N, l, B, H); // vector meson
	totw = hbar/tau0B;
	return pw/totw;
}


// decay BR: semileptonic case, 2 body

Double_t decayBR_2body_semilept(int idl, int idH, Double_t mN, Double_t tau0mN){
	
	Double_t BR;
	Double_t pw, totw;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	//Declare the HNL
	HNL N = HNL("HNL", mN, U2, mixes_with);
	N.setMajorana(majorana);
	
	Meson H; Lepton l;
	int i,j;
	
	switch(idH){
		case 111: H = pi; j = 0; break;
	}
		
	switch(idl){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
			
	pw = pw_neutral_pseudoscalar_mesons(cfg, l, H, N);
	totw = hbar/tau0mN;
	return pw/totw;
}








   

