# This script produces figure 2 from Bhattacharyya and Nityananda (2008).
# Use it to double check the PB found looks realistic.

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, tan, atan, sqrt, pi
from scipy.optimize import fsolve
import matplotlib as mpl

class trans(object):
  ecc = 0.0
  def f(self,E,M):
    return E-self.ecc*sin(E)-M
  def getE(self,M):
    return fsolve(self.f,0,args=(M))

def mass_func(pb_d,asini_c):
  msun  = 2.e30
  c     = 3.e8
  G     = 6.67e-11
  pb    = pb_d*24*60*60 
  asini = c*asini_c

  f = (4*pi**2/G)*(asini)**3/pb**2
  
  return f/msun

# input_params is a list of [p0,pb,a1,t0,om,ecc]
def per_vs_phi(input_params,sorted_phase,sorted_pers,sorted_errs,outfile="orbit_model.pdf"):

  p0,pb,a,t0,w,e = input_params
  p0_plot = f"{p0:.4f}"

  for sph,spe in zip(sorted_phase, sorted_pers):
    print(sph, spe)
  print(np.min(sorted_pers)-float(p0_plot), np.max(sorted_pers)-float(p0_plot))

  plt.errorbar(sorted_phase,np.array(sorted_pers)-float(p0_plot),yerr=sorted_errs,color='black',fmt='.',capsize=3,label='Measured Spin Period')

  plt.xlabel('Orbital Phase')
  plt.ylabel(f"Measured Spin Period $-$ {p0_plot} (ms)")

  phase, P_obs = orb_guess(p0,pb,e,w,a)

  if np.max(P_obs) > np.max(sorted_pers): ymax = np.max(P_obs)-float(p0_plot)
  else: ymax = np.max(sorted_pers)-float(p0_plot)

  if np.min(P_obs) < np.min(sorted_pers): ymin = np.min(P_obs)-float(p0_plot)
  else: ymin = np.min(sorted_pers)-float(p0_plot)

  y_values = (np.array(P_obs)-float(p0_plot))
  plt.plot(phase/(2*pi),np.array(P_obs)-float(p0_plot),ls='--',color='gray',label='Binary Model')
  plt.xlim([0,1])
  plt.ylim([1.1*ymin,1.1*ymax])

  f = mass_func(pb,a)
  print("")
  print("MASS FUNCTION: ",f)
  plt.legend()
  plt.savefig(outfile,format='pdf',bbox_inches='tight')


def orb_guess(p_0,pb_day,ecc,w_deg,asini_lts):

  func = trans()
  func.ecc = ecc
  c = 3.0e8

  phase = np.linspace(0,1.0,num=1000)*2*pi       # This is equivalent to "mean anomaly" (M) in HBOPA
  w = w_deg*pi/180.0
  pb = pb_day*24*60*60
  E   = []
  A_T = []
  V_1 = []
  P_obs = []
  angfreq_orb = 2.0*pi/pb
  inc = pi/2.0

  for i in range(len(phase)):
    E.append(func.getE(phase[i]))
    A_T.append(2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(E[i]/2.0)))       # atan returns A_T in rad
    V_1.append(angfreq_orb*(c*asini_lts*sin(inc)/sqrt(1.0-ecc**2))*(cos(w+A_T[i])+ecc*cos(w)))
    P_obs.append(p_0*(1.0+V_1[i]/c))

  return phase,P_obs

#### MAIN ####

# Set some plotting parameters
plt.rc("font",**{"family":"serif","serif":["Computer Modern Roman"]})
fig_width = 6.0
fig_height = fig_width*0.8
fig_size = [fig_width,fig_height]
params = {"backend": "pdf",
          "font.size"      : 14,
          "axes.labelsize" : 14,
          "legend.fontsize": 12,
          "xtick.labelsize": 14,
          "ytick.labelsize": 14,
          "text.usetex"    : True,
          "figure.figsize" : fig_size,
          "axes.unicode_minus": True}
mpl.rcParams.update(params)

# earlier version had all sizes modified by *fig_width/8.5

# Hard-coded for 1017
fname = "data/1017-mjd-p0-err.txt"
ecc = 0.227750
om = 60.012
t0 = 57545.748703
pb = 8.98397
a1 = 26.15660
p0 = 83.15253
input_params = [p0,pb,a1,t0,om,ecc]

# Read data file.
data = np.loadtxt(fname,dtype='float')
epoc = data[:,0]
pers = data[:,1] #*1000.0
errs = data[:,2] #*1000.0

# Convert epochs to orbital "phase" (mean anomaly)
phase = ((epoc-t0)/pb)%1.0
inds = np.argsort(phase)
sort_pers = pers[inds]
sort_phase = phase[inds]
sort_errs = errs[inds]
per_vs_phi(input_params,sort_phase,sort_pers,sort_errs,outfile="1017_orb_model.pdf")
