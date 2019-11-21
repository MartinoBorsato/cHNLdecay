from subprocess import Popen, PIPE

#def get_prod_BR(HNLmass_GeV, HNLlifetime_ns, B_ID, meson_ID, lepton_ID, cHNLdecayDIR='../cHNLdecay/', verbose=False):
    #HNLmass_MeV = HNLmass_GeV*1e3
    #output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--daughterMesonID', str(meson_ID), '--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
    #out = output.stdout.read()
    #out_BR = float(out.split()[0])
    #out_U2 = float(out.split()[1])
    #if verbose:
        #print("\n------\nINFO: BR calculator output: ", out_BR, "\n-------\n")

    #return out_BR, out_U2

#def get_prod_BR_U2(HNLmass_GeV, U2, B_ID, meson_ID, lepton_ID, cHNLdecayDIR='../cHNLdecay/', verbose=False):
    #HNLmass_MeV = HNLmass_GeV*1e3
    #output = Popen(['./cHNLdecay', '--mainmode', '2', '--BmesonID', str(B_ID), '--daughterMesonID', str(meson_ID), '--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
    #out = output.stdout.read()
    #out_BR = float(out.split()[0])
    #out_tau = float(out.split()[1])
    #if verbose:
        #print("\n------\nINFO: BR calculator output: ", out_BR, "\n-------\n")

    #return out_BR, out_tau

def get_prod_BR(HNLmass_GeV, HNLlifetime_ns, B_ID, meson_ID, lepton_ID, cHNLdecayDIR='../cHNLdecay/'): 
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--PrimaryMesonID', str(meson_ID),\
					'--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");

	out_BR = float(out.split()[0]);
	return out_BR;
	
def get_decay_BR_lepton_meson(HNLmass_GeV, U2, lepton_ID, meson_ID, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '2', '--SecondaryMesonID', str(meson_ID), '--LeptonA_ID', str(lepton_ID),\
					'--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_BR = float(out.split()[0]);
	
	return out_BR;
	#return 0;
	
def get_decay_BR_lepton_lepton_neutrino(HNLmass_GeV, HNLlifetime_ns, leptonA_ID, leptonB_ID, neutrinoB_ID, generations, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '3', '--LeptonA_ID', str(leptonA_ID),\
					'--LeptonB_ID', str(leptonB_ID), '--NeutrinoB_ID', str(neutrinoB_ID),\
					'--generations', str(generations), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_BR = float(out.split()[0]);
	
	return out_BR;
	
	
	
def get_Umu2(HNLmass_GeV, HNLlifetime_ns, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '6', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_Umu2 = float(out.split()[0]);
	
	return out_Umu2;
	

def get_tau0ns(HNLmass_GeV, U2, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '7', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print "\n------\nINFO: BR calculator output: ", out, "\n-------\n";
	
	out_tau0ns = float(out.split()[0]);
	
	return out_tau0ns;


def get_GammaOverGammaInv(HNLmass_GeV, U2, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	output = Popen(['./cHNLdecay', '--mainmode', '8', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print "\n------\nINFO: BR calculator output: ", out, "\n-------\n";
	
	out_GoGinv = float(out.split()[0]);
	
	return out_GoGinv;

def get_GammaMesonOverGammaQuarks(HNLmass_GeV, U2, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	output = Popen(['./cHNLdecay', '--mainmode', '9', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print "\n------\nINFO: BR calculator output: ", out, "\n-------\n";
	
	out_Gratio = float(out.split()[0]);
	
	return out_Gratio;



import matplotlib
matplotlib.use('Agg')
#from matplotlib import pyplot as plt
import pylab as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys

def plot_prod_BR():
	
	BR11 = []; BR12 = []; BR13 = []; BR14 = []
	BR20 = []; BR21 = []; BR22 = []; BR23 = []; BR24 = []; 
	BR31 = []; BR32 = []; BR33 = []; BR34 = [];
	BR41 = [];
	M = (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)
	# Fixed lifetime 1ns
	lifetime=50 #1ns
	for mass in M:
		
		BR11.append(get_prod_BR(mass, lifetime, 511, 211, 13))
		BR12.append(get_prod_BR(mass, lifetime, 511, 213, 13))
		BR13.append(get_prod_BR(mass, lifetime, 511, 411, 13))
		BR14.append(get_prod_BR(mass, lifetime, 511, 413, 13))
		
		BR21.append(get_prod_BR(mass, lifetime, 521, 0, 13))
		BR22.append(get_prod_BR(mass, lifetime, 521, 111, 13))
		BR23.append(get_prod_BR(mass, lifetime, 521, 113, 13))
		BR24.append(get_prod_BR(mass, lifetime, 521, 421, 13))
		BR20.append(get_prod_BR(mass, lifetime, 521, 423, 13))
		
		BR31.append(get_prod_BR(mass, lifetime, 531, 321, 13))
		BR32.append(get_prod_BR(mass, lifetime, 531, 323, 13))
		BR33.append(get_prod_BR(mass, lifetime, 531, 431, 13))
		BR34.append(get_prod_BR(mass, lifetime, 531, 431, 13))
		
		BR41.append(get_prod_BR(mass, lifetime, 541, 0, 13))
		
		
		
		
		#np.append(BR3, 0.5, get_prod_BR(m, 1, 511, 0, 13))
	
	plt.figure()
	
	plt.yscale('log')
	plt.title(r'BR prod, for fixed $\tau_N$ = '+str(lifetime)+' ns')
	plt.xlabel("$m_N$ [GeV]", fontsize=14)
	plt.ylabel('BR', fontsize=14)
	
	
	plt.plot(M, BR11, color='sienna', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow \pi^- \mu^+ N$')
	plt.plot(M, BR12, color='orangered', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow \rho^- \mu^+ N$')
	plt.plot(M, BR13, color='tomato', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow D^- \mu^+ N$')
	plt.plot(M, BR14, color='darkred', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow D^{-*} \mu^+ N$')


	plt.plot(M, BR24, color='green', linestyle = 'dashed',linewidth = 0.8,label=r'$B^+ \rightarrow {\bar{D}^0}^* \mu^+ N$')
	plt.plot(M, BR23, color='blueviolet', linestyle = 'dashed',linewidth = 0.8,label=r'$B^+ \rightarrow {\bar{D}^0} \mu^+ N$')
	plt.plot(M, BR21, color='rebeccapurple', linestyle = 'dashed',linewidth = 0.8,label=r'$B^+ \rightarrow \pi^0 \mu^+ N$')
	plt.plot(M, BR22, color='mediumpurple', linestyle = 'dashed',linewidth = 0.8,label=r'$B^+ \rightarrow \rho^0 \mu^+ N$')

	plt.plot(M, BR20, color='slateblue', linestyle = 'solid',linewidth = 0.8,label=r'$B^+ \rightarrow \mu^+ N$')


	plt.plot(M, BR33, color='steelblue', linestyle = 'dashed', linewidth = 0.8,label=r'$B_s^0 \rightarrow D_s^- \mu^+ N$')
	plt.plot(M, BR34, color='aqua', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0_s \rightarrow {D_s^-}^* \mu^+ N$')
	plt.plot(M, BR31, color='teal', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0_s \rightarrow K^- \mu^+ N$')

	plt.plot(M, BR41, color='teal', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0_s \rightarrow K^- \mu^+ N$')

	
	plt.legend()
	plt.savefig('prod_BRs_check.pdf');

#def plot_decay_BR():
	
	#BR1 = [];
	#BR2 = [];
	#BR3 = [];
	#M = (0.5, 1, 2, 4, 6)
	## Fixed lifetime 1ns
	#lifetime=10 #1ns
	#for m in M:
		#BR1.append(get_decay_BR_lepton_meson(m, lifetime, 13, 211))
		##BR2.append(get_decay_BR_lepton_meson(m, lifetime, 13, 411))
		##np.append(BR3, 0.5, get_prod_BR(m, 1, 511, 0, 13))
	
	#plt.figure()
	#plt.yscale('log')
	#plt.title(r'BR prod, for fixed $\tau_N$ = '+str(lifetime)+' ns')
	#plt.xlabel("$m_N$ [GeV]", fontsize=14)
	#plt.ylabel('BR', fontsize=14)
	#plt.plot(M, BR1, color='sienna', linestyle = 'dashed',linewidth = 0.8,label=r'$N \rightarrow \mu^\pm \pi^\mp $')
	##plt.plot(M, BR2, color='orangered', linestyle = 'dashed',linewidth = 0.8,label=r'$N \rightarrow \mu^\pm K^\mp $')
	#plt.legend()
	#plt.savefig('decay_BRs_check.pdf');
	
#def plot_lifetime():
	
	#U2 = 1;
	#mN = np.linspace(0.5, 6.5, 20); #GeV
	#y=[]
	#for i in mN:
		#print i;
		#y.append(get_tau0ns(i,U2)*1e-9) #ns to s

	
	#plt.figure()
	#plt.xscale('log')
	#plt.yscale('log')
	#plt.title(r'lifetime, $U_{\mu N}^2$ = '+str(U2))
	#plt.xlabel("$m_N$ [GeV]", fontsize=14)
	#plt.ylabel('$tau_{0}$ [s]', fontsize=14)
	#plt.plot(mN, y, color='sienna', linestyle = 'dashed',linewidth = 0.8,label=r'')
	##plt.plot(M, BR2, color='orangered', linestyle = 'dashed',linewidth = 0.8,label=r'$N \rightarrow \mu^\pm K^\mp $')
	#plt.legend()
	#plt.savefig('lifetime_check.pdf');
	
