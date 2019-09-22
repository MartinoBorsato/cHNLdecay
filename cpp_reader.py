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

def get_prod_BR(HNLmass_GeV, HNLlifetime_ns, B_ID, meson_ID, lepton_ID, cHNLdecayDIR='../cHNLdecay/'): #FIXME need to add cwd
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--PrimaryMesonID', str(meson_ID),\
					'--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");

	out_BR = float(out.split()[0]);
	return out_BR;
	
def get_decay_BR_lepton_meson(HNLmass_GeV, HNLlifetime_ns, lepton_ID, meson_ID, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '2', '--SecondaryMesonID', str(meson_ID), '--LeptonA_ID', str(lepton_ID),\
					'--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_BR = float(out.split()[0]);
	
	return out_BR;
	#return 0;
	
def get_decay_BR_lepton_lepton_neutrino(HNLmass_GeV, HNLlifetime_ns, leptonA_ID, leptonB_ID, neutrinoB_ID, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '3', '--LeptonA_ID', str(leptonA_ID),\
					'--LeptonB_ID', str(leptonB_ID), '--NeutrinoB_ID', str(neutrinoB_ID),\
					'--generations', str(leptonA_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_BR = float(out.split()[0]);
	
	return out_BR;
	
	
	
def get_Umu2(HNLmass_GeV, HNLlifetime_ns, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '4', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_Umu2 = float(out.split()[0]);
	
	return out_Umu2;

def get_lifetime(HNLmass_GeV, U2, cHNLdecayDIR='../cHNLdecay/'):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '5', \
					'--generations','13', '--mass', str(HNLmass_MeV), '--angle', str(U2)], stdout=PIPE, cwd=cHNLdecayDIR)
	out = output.stdout.read()
	#print("\n------\nINFO: BR calculator output: ", out, "\n-------\n");
	
	out_lifetime = float(out.split()[0]);
	
	return out_lifetime;

