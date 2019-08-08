from subprocess import Popen, PIPE


def get_prod_BR(HNLmass_GeV, HNLlifetime_ns, B_ID, meson_ID, lepton_ID):
	
	HNLmass_MeV = HNLmass_GeV*1e3;
	
	output = Popen(['./cHNLdecay', '--mainmode', '1', '--BmesonID', str(B_ID), '--daughterMesonID', str(meson_ID), '--generations', str(lepton_ID), '--mass', str(HNLmass_MeV), '--lifetime-ns', str(HNLlifetime_ns)], stdout=PIPE)
	out = output.stdout.read()
	out_BR = float(out.split()[0]);
	print "\n------\nINFO: BR calculator output: ", out_BR, "\n-------\n";
	return out_BR;

import matplotlib
matplotlib.use('Agg')
#from matplotlib import pyplot as plt
import pylab as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys

def plot_prod_BR():
	
	BR1 = [];
	BR2 = [];
	BR3 = [];
	M = (0.5, 1, 2, 4)
	# Fixed lifetime 1ns
	lifetime=10 #1ns
	for m in M:
		BR1.append(get_prod_BR(m, lifetime, 521, 0, 13))
		BR2.append(get_prod_BR(m, lifetime, 521, 111, 13))
		#np.append(BR3, 0.5, get_prod_BR(m, 1, 511, 0, 13))
	
	plt.figure()
	plt.yscale('log')
	plt.title(r'BR prod, for fixed $\tau_N$ = '+str(lifetime)+' ns')
	plt.xlabel("$m_N$ [GeV]", fontsize=14)
	plt.ylabel('BR', fontsize=14)
	plt.plot(M, BR1, color='sienna', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow \mu^+ N$')
	plt.plot(M, BR2, color='orangered', linestyle = 'dashed',linewidth = 0.8,label=r'$B^0 \rightarrow \pi^- \mu^+ N$')
	plt.legend()
	plt.savefig('prod_BRs_check.pdf');
	
