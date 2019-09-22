
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

Double_t M_LAMBDAb = 5620.2;
Double_t M_PROTON  = 938.272;
Double_t M_LAMBDAc = 2286.46;
Double_t M_B = 5279.29;
Double_t M_Bc = 6275.1;
Double_t M_pi = 139.57018;

// Tables

Double_t get_fhhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.76;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.295;
	else if(hID == 531 and hpID == 433) f=0.95;
	else if(hID == 531 and hpID == 323) f=0.291;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_V" << std::endl;
	return f; 
}

Double_t get_fhhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.69;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.231;
	else if(hID == 531 and hpID == 433) f=0.67;
	else if(hID == 531 and hpID == 323) f=0.289;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A0" << std::endl;
	return f; 
}
Double_t get_fhhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.66;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.269;
	else if(hID == 531 and hpID == 433) f=0.70;
	else if(hID == 531 and hpID == 323) f=0.287;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A1" << std::endl;
	return f; 
}
Double_t get_fhhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.62;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.282;
	else if(hID == 531 and hpID == 433) f=0.75;
	else if(hID == 531 and hpID == 323) f=0.286;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A2" << std::endl;
	return f; 
}

Double_t get_sigmahhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.57;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.875;
	else if(hID == 531 and hpID == 433) sigma=0.372;
	else if(hID == 531 and hpID == 323) sigma=0.516; // changed
	else std::cerr << "ERROR: in get_sigmahhp_V" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.59;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.796;
	else if(hID == 531 and hpID == 433) sigma=0.350;
	else if(hID == 531 and hpID == 323) sigma=-0.383;
	else std::cerr << "ERROR: in get_sigmahhp_A0" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.78;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.54;
	else if(hID == 531 and hpID == 433) sigma=0.463;
	else if(hID == 531 and hpID == 323) sigma=0.;
	else std::cerr << "ERROR: in get_sigmahhp_A1" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=1.40;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=1.34;
	else if(hID == 531 and hpID == 433) sigma=1.04;
	else if(hID == 531 and hpID == 323) sigma=1.05;
	else std::cerr << "ERROR: in get_sigmahhp_A2" << std::endl;
	return sigma; 
}
//
Double_t get_xihhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.;
	else if(hID == 531 and hpID == 433) xi=0.561;
	else if(hID == 531 and hpID == 323) xi=2.1;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.055;
	else if(hID == 531 and hpID == 433) xi=0.600;
	else if(hID == 531 and hpID == 323) xi=1.58;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.;
	else if(hID == 531 and hpID == 433) xi=0.510;
	else if(hID == 531 and hpID == 323) xi=1.06;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.41;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=-0.21;
	else if(hID == 531 and hpID == 433) xi=0.070;
	else if(hID == 531 and hpID == 323) xi=-0.074;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}


Double_t get_MhP(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t M(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) M=6275.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) M=5279.;
	else if(hID == 531 and hpID == 433) M=6275.;
	else if(hID == 531 and hpID == 323) M=5367.;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_MhP" << std::endl;
	return M; 
}
Double_t get_MhV(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t M(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) M=6331.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) M=5325.;
	else if(hID == 531 and hpID == 433) M=6331.;
	else if(hID == 531 and hpID == 323) M=5415.;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_MhV" << std::endl;
	return M; 
}
 
 


// Functions

Double_t compute_FF_V(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_V(h, hp);
	Double_t sigmahhp = get_fhhp_V(h, hp);
	Double_t xihhp = get_xihhp_V(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-(q2/pow(MV,2)))*(1-sigmahhp*(q2/pow(MV,2)) - xihhp*pow(q2,2)/pow(MV,4)));
	//if((hID == 531 and hpID == 323))std::cout<<"V: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A0(Meson h, Meson hp, Double_t q2){
	
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A0(h, hp);
	Double_t sigmahhp = get_fhhp_A0(h, hp);
	Double_t xihhp = get_xihhp_A0(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-q2/pow(MP,2))*(1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A0: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A1(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A1(h, hp);
	Double_t sigmahhp = get_fhhp_A1(h, hp);
	Double_t xihhp = get_xihhp_A1(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A1: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A2(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A2(h, hp);
	Double_t sigmahhp = get_fhhp_A2(h, hp);
	Double_t xihhp = get_xihhp_A2(h, hp);
	Double_t MP = get_MhP( h, hp);
	Double_t MV = get_MhV( h, hp);
	
	res = (fhhP)/((1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A2: "<<res<<std::endl;
	return res;
}


// -------------
Double_t compute_FF_g(Meson h, Meson hp, Double_t q2){
	Double_t res, V;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	V = compute_FF_V(h, hp, q2);
	res = V/(mh+mhp);
	return res;
}

Double_t compute_FF_aplus(Meson h, Meson hp, Double_t q2){
	Double_t res, A2;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A2 = compute_FF_A2(h, hp, q2);
	res = -1*A2/(mh+mhp);
	return res;
}

Double_t compute_FF_f(Meson h, Meson hp, Double_t q2){
	Double_t res, A1;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A1 = compute_FF_A1(h, hp, q2);
	res = (mh+mhp)*A1;
	return res;
}

Double_t compute_FF_aminus(Meson h, Meson hp, Double_t q2){
	Double_t res, A0;
	Double_t aplus = compute_FF_aplus(h, hp, q2);
	Double_t f = compute_FF_f(h, hp, q2);
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A0 = compute_FF_A0(h, hp, q2);
	res = (2*mhp*A0-f-(pow(mh,2)-pow(mhp,2))*aplus)/q2;
	return res;
}



Double_t compute_ffactor(Meson h, Meson hp, Double_t q2, int opt) {
	
	// Uses equation (C.16) of arXiv paper
	
	Double_t 	Mpole, mh, mhp, a0, a1, a2, z, tplus, t0, factor, sum;
	Int_t 		IDh, IDhp;
	
	mh =  h.getMass(); 		mhp =  hp.getMass();
	IDh = h.getPdgId();		IDhp = hp.getPdgId();
	
	switch(opt){
		
		case 1:				//form factor fplus
				
			switch(IDhp){	
				case 411:			//D
				//case 413:			//Dstar
				case 421:			//D0bar
				//case 423:			//D0starbar
				case 431:			//Ds
				//case 433:			//Dsstar
					//Mpole = inf;
					a0 = 0.909;
					a1 = -7.11;
					a2 = 66;
					factor = 1.;
					break;
				case 321:			//K
				//case 323:			//Kstar
					Mpole = 5325;
					a0 = 0.360;
					a1 = -0.828;
					a2 = 1.1;
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
				case 211:			//pi
				case 111:			//pi0
					Mpole = 5325;
					a0 = 0.404;
					a1 = -0.68;
					a2 = -0.86;	
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
			}
			break;
		case 0: // form factor f0
			switch(IDhp){	
				case 411:			//D
				case 413:			//Dstar
				case 421:			//D0bar
				case 423:			//D0starbar
				case 431:			//Ds
				//  case 433:			//Dsstar
					//Mpole = inf;
					a0 = 0.794;
					a1 = -2.45;
					a2 = 33;
					factor = 1.;
					break;
				case 321:			//K
				//case 323:			//Kstar
					Mpole = 5650;
					a0 = 0.233;
					a1 = 0.197;
					a2 = 0.18;
					factor = 1/(1-q2/pow(Mpole, 2));
					//std::cout<<"factor " <<std::endl;
					break;
				case 211:			//pi
				case 111:			//pi0
					Mpole = 5650;
					a0 = 0.490;
					a1 = -1.61;
					a2 = 0.93;	
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
					//std::cout<<"factor " <<std::endl;
			}
			break;
				
		case 2:		// g
			//std::cout<<"g"<<compute_FF_g(h, hp, q2)<<std::endl;
			return compute_FF_g(h, hp, q2);
			//return 1.;
			break;
		case 3:		// f	
			//std::cout<<"f:"<<compute_FF_f(h, hp, q2)<<std::endl;
			//return 1.;
			return compute_FF_f(h, hp, q2);
			break;
		case 4:		// aplus	
			//std::cout<<"aplus:"<<compute_FF_aplus(h, hp, q2)<<std::endl;
			//return 1.;
			return compute_FF_aplus(h, hp, q2);
			break;
		case 5:		// aminus	
			//return 1.;
			//std::cout<<"aminus:"<<compute_FF_aminus(h, hp, q2)<<std::endl;
			return compute_FF_aminus(h, hp, q2);
			break;
	
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

// Form factors hyperons

std::vector<Double_t> acoeffs_Lambdab_proton() // table X of [1503.01421]
{
	return { 0.42512661779076, -0.70879824786202, 0.89250038017833, 
			 0.41441295090913, -1.041709976979, 1.9259425862863, 
			 0.52144673308125, -0.82467696338466, 0.76090923603286,
			 0.38886025531083,  -1.0725956052846, 1.9857217253182,
			 0.44189326677855,  -0.86494888515889, 0.99692910373979,
			 0.38886025531083, -1.0839467598681, 1.4521580817004 };

	/* { a0_fplus, 	a1_fplus,	a2_fplus  
		 a0_f0,		a1_f0, 		a2_f0  
		 a0_fperp,  a1_fperp, 	a2_fperp 
		 a1_gplus,	a2_gplus 
		a0_g0, 		a1_g0, 		a2_g0 
		a1_gperp, a2_gperp

* */
}


std::vector<Double_t> acoeffs_Lambdab_Lambdac() // table X of [1503.01421]
{
	return { 0.81034719382496, -4.7475223120166, 0.78618045366394, // fplus
			 0.73889965204015, -4.5633603364873, 2.7050827605253, //f0
			 1.0939236137642, -6.4413374884146, 2.3161239651561, //fperp
			 0.68483287500112, -4.3788469218438, 1.2813964795243, //gplus
			 0.74078109306203, -4.3860778472363, 1.3375076880839, //g0
			 0.68483287500112, -4.6270436667368, 1.6143367954872  //gperp
			};
}


std::vector<Double_t> Mpole_Lambdab_proton(){
	return{ M_B + 46,  M_B + 377,  M_B + 46, M_B + 427, M_B, M_B + 427 };	
	// fplus, f0, fperp, gplus, g0, gplus
}
std::vector<Double_t> Mpole_Lambdab_Lambdac(){
	return{ M_Bc + 56,  M_Bc + 449,  M_Bc + 56, M_Bc + 492, M_Bc, M_Bc + 492 };	
	// fplus, f0, fperp, gplus, g0, gplus
}
std::vector<Double_t> tplus_Lambdab_proton(){
	return{ pow(M_B+M_pi,2), pow(M_B+M_pi,2), pow(M_B+M_pi,2),
			pow(M_B+M_pi,2), pow(M_B+M_pi,2), pow(M_B+M_pi,2)
			};
	// fplus, f0, fperp, gplus, g0, gplus
}
std::vector<Double_t> tplus_Lambdab_Lambdac(){
	std::vector<Double_t> result;
	for(size_t i(0); i<6; ++i) {
		result.push_back(Mpole_Lambdab_Lambdac()[i]);
	}
	return result;
	// fplus, f0, fperp, gplus, g0, gplus
}


Double_t t0_proton(){
	return pow(M_LAMBDAb - M_PROTON,2);
}
Double_t t0_Lambdac(){
	return pow(M_LAMBDAb - M_LAMBDAc,2);
}
Double_t zfunction_Lambdab_proton(Double_t tplus, Double_t q2){
	return (sqrt(tplus - q2) - sqrt(tplus - t0_proton()))/(sqrt(tplus-q2)+sqrt(tplus -t0_proton()));
}
Double_t zfunction_Lambdab_Lambdac(Double_t tplus, Double_t q2){
	return (sqrt(tplus - q2) - sqrt(tplus - t0_Lambdac()))/(sqrt(tplus-q2)+sqrt(tplus -t0_Lambdac()));
}

std::vector<Double_t> z_Lambdab_proton(Double_t q2){
	
	std::vector<Double_t> result;
	Double_t tplus; 
	for(size_t i(0); i<6; ++i){
		tplus = tplus_Lambdab_proton()[i];
		result.push_back(zfunction_Lambdab_proton(tplus, q2));
		//std::cout<<"z["<<i<<"]"<<zfunction_Lambdab_Lambdac(tplus, q2)<<std::endl;
	}
	
	return result;
	// fplus, f0, fperp, gplus, g0, gplus
}

std::vector<Double_t> z_Lambdab_Lambdac(Double_t q2){
	
	std::vector<Double_t> result;
	Double_t tplus; 
	for(size_t i(0); i<6; ++i){
		tplus = tplus_Lambdab_Lambdac()[i];
		result.push_back(zfunction_Lambdab_Lambdac(tplus, q2));
		//std::cout<<"z["<<i<<"]"<<zfunction_Lambdab_Lambdac(tplus, q2)<<std::endl;
	}
	
	return result;
	// fplus, f0, fperp, gplus, g0, gplus
}


//// eq. (81)

Double_t f_Lambdab_proton(Double_t q2, int idx){
	
	Double_t Mpole, mh, mhp, a0, a1, a2, z, tplus, t0, factor, sum;
	
	
	//mh = M_LAMBDAb;
	//mhp = M_PROTON;
	
	
	
	//t0_proton		= pow(mh - mhp, 2); 
	tplus 	= tplus_Lambdab_proton()[idx];
	//std::cout<<"tplus" << tplus << std::endl;
	Mpole 	= Mpole_Lambdab_proton()[idx];
	//std::cout<<"Mpole" << Mpole << std::endl;
	z 		= z_Lambdab_proton(q2)[idx];
	//std::cout<<"z" << z << std::endl;
	a0		= acoeffs_Lambdab_proton()[idx*3];
	//std::cout<<"a0[" << idx*3 << "] = " << a0 << std::endl;
	a1		= acoeffs_Lambdab_proton()[idx*3+1];
	//std::cout<<"a1" << a1 << std::endl;
	//std::cout<<"a1[" << idx*3+1 << "] = " << a1 << std::endl;
	a2		= acoeffs_Lambdab_proton()[idx*3+2];
	//std::cout<<"a2[" << idx*3+2 << "] = " << a2 << std::endl;
	
	factor = 1/(1-q2/pow(Mpole, 2));
	//std::cout<<"factor" << factor << std::endl;
	
	sum  = a0;
	sum += a1*z;
	sum += a2*pow(z,2);
	
	//std::cout<<"sum" << sum << std::endl;
	
	//std::cout<<"factor*sum" << factor*sum << std::endl;
	return factor*sum;
	
}

Double_t f_Lambdab_Lambdac(Double_t q2, int idx){
	
	Double_t Mpole, mh, mhp, a0, a1, a2, z, tplus, t0, factor, sum;
	
	
	mh = M_LAMBDAb;
	mhp = M_LAMBDAc;
	
	
	
	//t0		= pow(mh - mhp, 2); 
	tplus 	= tplus_Lambdab_Lambdac()[idx];
	Mpole 	= Mpole_Lambdab_Lambdac()[idx];
	z 		= z_Lambdab_Lambdac(q2)[idx];
	std::cout<<"z" << z << std::endl;
	a0		= acoeffs_Lambdab_Lambdac()[idx*3];
	a1		= acoeffs_Lambdab_Lambdac()[idx*3+1];
	a2		= acoeffs_Lambdab_Lambdac()[idx*3+2];
	
	factor = 1/(1-q2/pow(Mpole, 2));
	
	sum  = a0;
	sum += a1*z;
	sum += a2*pow(z,2);
	
	return factor*sum;
}

////



Double_t f2V_Lambdab_proton(Double_t q2){ 
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	
	Double_t CplusV = q2/(mh*(mh+mhp));
	Double_t CperpV = (mh+mhp)/mh;
	/*
	std::cout << "f_Lambdab_proton(q2,0)" << f_Lambdab_proton(q2,0) << std::endl;
	std::cout << "f_Lambdab_proton(q2,1)" << f_Lambdab_proton(q2,1) << std::endl;
	std::cout << "f_Lambdab_proton(q2,2)" << f_Lambdab_proton(q2,2) << std::endl;
	std::cout << "f_Lambdab_proton(q2,3)" << f_Lambdab_proton(q2,3) << std::endl;
	std::cout << "f_Lambdab_proton(q2,4)" << f_Lambdab_proton(q2,4) << std::endl;
	std::cout << "f_Lambdab_proton(q2,5)" << f_Lambdab_proton(q2,5) << std::endl;
	*/
	
	//return (mh*(mh+mhp))/(q2-pow(mh+mhp,2)) * (f_Lambdab_proton(q2,0)-f_Lambdab_proton(q2,2));
	return (CplusV-CperpV) * (f_Lambdab_proton(q2,0)-f_Lambdab_proton(q2,2));
}
Double_t f1V_Lambdab_proton(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	
	Double_t CplusV = q2/(mh*(mh+mhp));
	
	return f_Lambdab_proton(q2,0) - CplusV*f2V_Lambdab_proton(q2);
}
Double_t f3V_Lambdab_proton(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	
	Double_t C0V = q2/(mh*(mh-mhp));
	//return (f_Lambdab_proton(q2,1) - f1V_Lambdab_proton(q2))*(mh*(mh-mhp))/q2; 
	return (f_Lambdab_proton(q2,1) - f1V_Lambdab_proton(q2))/C0V; 
}

Double_t f2A_Lambdab_proton(Double_t q2){ 
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	//std::cout << "q2 = " << q2 << std::endl;
	//std::cout<<"(mh*(mh-mhp))/(-q2+pow(mh-mhp,2))" << (mh*(mh-mhp))/(-q2+pow(mh-mhp,2)) << std::endl;
	
	Double_t CplusA = q2/(mh*(mh-mhp));
	Double_t CperpA = (mh-mhp)/mh;
	
	//return (mh*(mh-mhp))/(-q2+pow(mh-mhp,2)) * (f_Lambdab_proton(q2,3)-f_Lambdab_proton(q2,5));
	return (-CplusA+CperpA) * (f_Lambdab_proton(q2,3)-f_Lambdab_proton(q2,5));
	
}
Double_t f1A_Lambdab_proton(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	Double_t CplusA = q2/(mh*(mh-mhp));
	
	//return f_Lambdab_proton(q2,0) - (-q2/(mh*(mh-mhp)))*f2A_Lambdab_proton(q2);
	return f_Lambdab_proton(q2,0) - CplusA*f2A_Lambdab_proton(q2);
}
Double_t f3A_Lambdab_proton(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_PROTON;
	Double_t C0V = q2/(mh*(mh+mhp));
	
	//return (f_Lambdab_proton(q2,4) - f1V_Lambdab_proton(q2))*(mh*(mh+mhp))/(-q2); 
	return (f_Lambdab_proton(q2,4) - f1V_Lambdab_proton(q2))/C0V; 

}

///////////////////////////////////////////////////////////////////////

Double_t f2V_Lambdab_Lambdac(Double_t q2){ // true = V, false = A
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	return (mh*(mh+mhp))/(q2-pow(mh+mhp,2)) * (f_Lambdab_Lambdac(q2,0)-f_Lambdab_Lambdac(q2,2));
}
Double_t f1V_Lambdab_Lambdac(Double_t q2){ // true = V, false = A
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	
	return f_Lambdab_Lambdac(q2,0) - (q2/(mh*(mh+mhp)))*f2V_Lambdab_Lambdac(q2);
}
Double_t f3V_Lambdab_Lambdac(Double_t q2){ // true = V, false = A
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	
	return (f_Lambdab_Lambdac(q2,1) - f1V_Lambdab_Lambdac(q2))*(mh*(mh-mhp))/q2; 
}
Double_t f2A_Lambdab_Lambdac(Double_t q2){ 
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	return (mh*(mh-mhp))/(-q2+pow(mh-mhp,2)) * (f_Lambdab_Lambdac(q2,3)-f_Lambdab_Lambdac(q2,5));
}
Double_t f1A_Lambdab_Lambdac(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	
	return f_Lambdab_Lambdac(q2,0) - (-q2/(mh*(mh-mhp)))*f2A_Lambdab_Lambdac(q2);
}
Double_t f3A_Lambdab_Lambdac(Double_t q2){ 
	
	Double_t mh = M_LAMBDAb;
	Double_t mhp = M_LAMBDAc;
	
	return (f_Lambdab_Lambdac(q2,4) - f1V_Lambdab_Lambdac(q2))*(mh*(mh+mhp))/(-q2); 
}


Double_t Sigma_hyp(int sgn, Double_t mh, Double_t mhp){
	return pow(mh, 2) +sgn*pow(mhp, 2);
}
Double_t Delta_hyp(int sgn, Double_t mh, Double_t mhp){
	return mh + sgn*mhp; 
}

/////////////////

// Equations 

Double_t a1(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V){
	Double_t Sigmap = Sigma_hyp(+1, mh, mhp); Double_t Sigmam = Sigma_hyp(-1, mh, mhp);
	Double_t Deltap = Delta_hyp(+1, mh, mhp); Double_t Deltam = Delta_hyp(-1, mh, mhp);
	
	/*
	std::cout << "Sigmap: " << Sigmap << std::endl;
	std::cout << "Deltap: " << Deltap << std::endl;
	std::cout << "Sigmam: " << Sigmam << std::endl;
	std::cout << "Deltam: " << Deltam << std::endl;
	*/
	Double_t Sigma_VA, Delta_VA; int sgn;
	if(V){
		sgn = +1;
		Sigma_VA = Sigmap;
		Delta_VA = Deltam;
	}
	else{
		sgn = -1;
		Sigma_VA = Sigmam;
		Delta_VA = Deltap;
	}
	
	/*
	std::cout << "------ a1 = " << pow(ml,2)*(q2*(pow(Sigmam,2) - 2*pow(mN,2)*Sigmap)
					 - 2*pow(q2,2)*(pow(mN,2) - 1*sgn*mh*mhp + pow(Delta_VA,2))
					 + 4*pow(mN,2)*pow(Sigmam, 2) 
					 + pow(q2, 3)
					 )
		 +(q2 - pow(mN,2))* (pow(mN,2)*(q2*Sigmap-2*Sigmam+pow(q2,2))
							 - q2*((pow(Delta_VA,2)-q2)*(Sigma_VA+2*q2))
							) << std::endl;
	*/
	return pow(ml,2)*(q2*(pow(Sigmam,2) - 2*pow(mN,2)*Sigmap)
					 - 2*pow(q2,2)*(pow(mN,2) - sgn*mh*mhp + pow(Delta_VA,2))
					 + 4*pow(mN,2)*pow(Sigmam, 2) 
					 + pow(q2, 3)
					 )
		 +(q2 - pow(mN,2))* (pow(mN,2)*(q2*Sigmap-2*Sigmam+pow(q2,2))
							 - q2*((pow(Delta_VA,2)-q2)*(pow(Sigma_VA,2)+2*pow(q2,2))) //pow(q2,2) instead of q2
							);
}

Double_t a2(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V){
	Double_t Sigmap = Sigma_hyp(+1, mh, mhp); Double_t Sigmam = Sigma_hyp(-1, mh, mhp);
	Double_t Deltap = Delta_hyp(+1, mh, mhp); Double_t Deltam = Delta_hyp(-1, mh, mhp);
	
	Double_t Sigma_VA, Delta_VA; int sgn;
	if(V){
		sgn=+1;
		Sigma_VA = Sigmap;
		Delta_VA = Deltam;
	}
	else{
		sgn=-1;
		Sigma_VA = Sigmam;
		Delta_VA = Deltap;
	}
	
	
	std::cout << "---- a2 part 1 :" << ( q2*(pow(ml,2) + pow(mN, 2)) 
			  +pow(pow(ml,2)-pow(mN,2),2)
			  -2*pow(q2,2) 
			  ) << std::endl;
	std::cout << "---- a2 part 2 :" << q2*( pow(Sigma_VA,2)+sgn*4*mh*mhp )
			  -2*(pow(Sigmam,2)+pow(q2,2)) << std::endl;		  
			  
	
	return  ( q2*(pow(ml,2) + pow(mN, 2)) 
			  +pow((pow(ml,2)-pow(mN,2)),2)
			  -2*pow(q2,2) 
			  )
			*(q2*( pow(Sigma_VA,2)+sgn*4*mh*mhp )
			  -2*(pow(Sigmam,2)+pow(q2,2))
			  );
}

Double_t a3(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V){
	Double_t Sigmap = Sigma_hyp(+1, mh, mhp); Double_t Sigmam = Sigma_hyp(-1, mh, mhp);
	Double_t Deltap = Delta_hyp(+1, mh, mhp); Double_t Deltam = Delta_hyp(-1, mh, mhp);
	
	Double_t Sigma_VA, Delta_VA; int sgn;
	if(V){
		sgn=+1;
		//Sigma_VA = Sigmap;
		Delta_VA = Deltap;
	}
	else{
		sgn=-1;
		//Sigma_VA = Sigmam;
		Delta_VA = Deltam;
	}
	return (pow(Delta_VA,2) - pow(q2,2)) * 
				 ( pow(ml,2)* (  q2 + 2*pow(mN,2) - pow(ml,2))
				  +pow(mN,2)*(q2 - pow(mN,2)) );
}
Double_t a12(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V){
	Double_t Sigmap = Sigma_hyp(+1, mh, mhp); Double_t Sigmam = Sigma_hyp(-1, mh, mhp);
	Double_t Deltap = Delta_hyp(+1, mh, mhp); Double_t Deltam = Delta_hyp(-1, mh, mhp);
	
	//Double_t Sigma_VA, Delta_VA; 
	if(V){
		return Sigmap*(Sigmam-q2)*
		( pow(ml,2)*(q2-2*pow(mN,2)) + pow(ml,4) + pow(mN,4) + pow(mN,2)*q2 - 2*pow(q2,2));
	}
	else{
		return Sigmam*(Sigmap-q2)*
		( pow(ml,2)*(q2-2*pow(mN,2)) + pow(ml,4) + pow(mN,4) + pow(mN,2)*q2 - 2*pow(q2,2));
	}	

}
Double_t a13(Double_t q2, Double_t mh, Double_t mhp, Double_t ml, Double_t mN, bool V){
	Double_t Sigmap = Sigma_hyp(+1, mh, mhp); Double_t Sigmam = Sigma_hyp(-1, mh, mhp);
	Double_t Deltap = Delta_hyp(+1, mh, mhp); Double_t Deltam = Delta_hyp(-1, mh, mhp);
	
	//Double_t Sigma_VA, Delta_VA; 
	if(V){
		return Sigmam*(Sigmap-q2)*
		( pow( pow(ml,2) - pow(mN,2) ,2) - q2*(pow(ml,2) + pow(mN,2)) );
	}
	else{
		return Sigmap*(Sigmam-q2)*
		( pow( pow(ml,2) - pow(mN,2) ,2) - q2*(pow(ml,2) + pow(mN,2)) );
	}
}


