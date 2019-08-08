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

double hbar=6.582119514e-21; //in MeV/s
std::vector<std::vector<Double_t>> VCKM{{0.97427, 0.22534, 0.00351}, 
										{0.22520, 0.97344, 0.04120}, 
										{0.00867, 0.04040, 0.999146}};

Double_t mN = 0;
Double_t U2 = 0;

bool majorana = true;

// variables used for searching for angle / parsing file for angles
bool isLoad = false;
//TString loadPath;
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


/** ERRORS MESSAGES **/

void test_value(Double_t var, Double_t minval, Double_t maxval, std::string varname){

	try {
		if(var>maxval or var<minval){
			std::string err_msg = "\n-------\nERROR: "+varname+" must outside allowed intervall ["+std::to_string(minval)+", "+std::to_string(maxval)+"]! \nPlease check your input mass[GeV] and lifetime[ns].\n-------\n";
			throw err_msg;
		}
	}
	catch (const std::string msg) {
		std::cerr << msg << std::endl;
	}
	return;
}

/** UTILS **/

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::vector<double> > file_to_vec(std::string filename){
	
	int nlines;		// n of points	
	int nrows(2); 	// x-y plot
	
	std::vector<std::vector<double> > tau1_vec(213,std::vector<double>(2, 0));
	std::ifstream file_tau1(filename);
	std::string line;
	int i1(0);
	while (std::getline(file_tau1, line)){

        std::vector<std::string> row_values;
        split(line, ' ', row_values);
        for (int i2(0); i2<row_values.size(); i2++){
			std::string v=row_values[i2];
			tau1_vec[i1][i2]=atof(v.c_str());
		}
		i1++;
	}
	return tau1_vec;
}

// Returns element closest to target in arr[] 


int getClosestIdx(Double_t val1, int idx1, Double_t val2, int idx2, Double_t target) 
{ 
    if (target - val1 >= val2 - target) 
        return idx2; 
    else
        return idx1; 
} 

int findClosestIdx(std::vector<Double_t> arr, Double_t target) 
{ 
	int n = arr.size();
    // Corner cases 
    if (target <= arr[0]) 
        //return arr[0]; 
        return 0;
    if (target >= arr[n - 1]) 
        //return arr[n - 1]; 
        return n-1;
  
    // Doing binary search 
    int i = 0, j = n, mid = 0; 
    while (i < j) { 
        mid = (i + j) / 2; 
  
        if (arr[mid] == target) 
            //return arr[mid]; 
            return mid;
  
        // If target is less than array element then search in left 
        if (target < arr[mid]) { 
  
            // If target is greater than previous to mid, return closest of two 
            if (mid > 0 && target > arr[mid - 1])
                return getClosestIdx(arr[mid - 1], mid-1, arr[mid], mid, target); 
  
            // Repeat for left half */
            j = mid; 
        } 
  
        // If target is greater than mid 
        else { 
            if (mid < n - 1 && target < arr[mid + 1]) 
                return getClosestIdx(arr[mid], mid, arr[mid + 1], mid+1, target); 
            // update i 
            i = mid + 1;  
        } 
    } 
    // Only single element left after search 
    //return arr[mid]; 
    return mid;
} 


/** FUNCTIONS **/

Double_t tau0_to_U2(Double_t mN, Double_t tau0mN){
	
	
	std::vector<std::vector<Double_t>> mNtau1_vec = file_to_vec("./HNLtau1_xsorted.csv");
	size_t s = mNtau1_vec.size();
	std::vector<Double_t> mN_vec(s,0);
	std::vector<Double_t> tau1_vec(s,0);
	
	for(size_t i(0); i<s; ++i){
		mN_vec[i] = mNtau1_vec[i][0];
		tau1_vec[i] = mNtau1_vec[i][1];
	}
	
	Double_t U2, tau1;
	int idxmN_;
	//computations
	
	//grab closer mN (convert in GeV to compare to the paper)
	idxmN_ = findClosestIdx(mN_vec, mN*1e-3); 	//mN in GeV
	tau1 = tau1_vec[idxmN_]*1e9; 				// tau in ns
	
	//proportionnality
	U2=tau1/tau0mN;
	return U2;
}


// production BR: leptonic case

Double_t prodBR_lept(int idB, int idl, Double_t mN, Double_t tau0mN){
	
	Double_t BR;
	Double_t pw, totw, tau0B;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	test_value(U2, 0., 1., "coupling U_{muN}^2");
	
	//std::cout << "U2(" << tau0mN << ") = " << U2 << std::endl;
	
	//Declare the HNL
	HNL N = HNL("HNL", mN, U2, mixes_with);
	N.setMajorana(majorana);
	
	Meson B; Lepton l;
	int i,j;
	
	switch(idB){
		case 521: B = Bp; i = 2; j = 0; tau0B = tau0Bp; break;
		case 541: B = Bc; i = 2; j = 1; tau0B = tau0Bc; break;
		default: std::cerr<<"ERROR: Meson ID not among the pre-programmed list of leptonic decays (521, 541)!"<<std::endl; return 1.;	
	}
		
	switch(idB){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
			
	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_leptonic(cfg, N, l, B);
	totw = hbar/tau0B;
	test_value(pw/totw, 0., 1., "Production branching ratio");
	return pw/totw;
}

// production BR: semileptonic case in pseudoscalar meson 

Double_t prodBR_semilept(int idB, int idl, int idH, Double_t mN, Double_t tau0mN){
	
	
	Double_t BR;
	Double_t pw, totw, tau0B;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	test_value(U2, 0., 1., "coupling U_{muN}^2");

	
	
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
		default: std::cerr<<"ERROR: Beauty meson ID not among the pre-programmed list of semileptonic decays (511, 521, 531, 541)!"<<std::endl;	
	}
	
	// try to do it with jkey
	switch(idH){
		case 111: H = pi; j = 0; break;
		default: std::cerr<<"ERROR: Meson ID not among the pre-programmed list (111)!"<<std::endl; return 1.;

	}
		
	switch(idB){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
	
	if(isP)	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_semileptonic(cfg, N, l, B, H); // pseudoscalar meson
	else 	pw = pow(VCKM[i][j],2)*pw_prodFromBmeson_semileptonic(cfg, N, l, B, H); // vector meson
	totw = hbar/tau0B;
	
	test_value(pw/totw, 0., 1., "Production branching ratio");
	
	return pw/totw;
}


// decay BR: semileptonic case, 2 body

Double_t decayBR_2body_semilept(int idl, int idH, Double_t mN, Double_t tau0mN){
	
	Double_t BR;
	Double_t pw, totw;
	Double_t U2 = tau0_to_U2(mN, tau0mN);
	
	test_value(U2, 0., 1., "coupling U_{muN}^2");
	
	//Declare the HNL
	HNL N = HNL("HNL", mN, U2, mixes_with);
	N.setMajorana(majorana);
	
	Meson H; Lepton l;
	int i,j;
	
	switch(idH){
		case 111: H = pi; j = 0; break;
		default: std::cerr<<"ERROR: Meson ID not among the pre-programmed list (111)!"<<std::endl;	
	}
		
	switch(idl){
		case 11: l = el; break;
		case 13: l = mu; break;
		case 15: l = tau; break;
	}
			
	pw = pw_neutral_pseudoscalar_mesons(cfg, l, H, N);
	totw = hbar/tau0mN;
	
	test_value(pw/totw, 0., 1., "Production branching ratio");
	return pw/totw;
}









   

