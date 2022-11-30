import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

data = np.loadtxt("DMX_vals_errs.txt", dtype=str, usecols=[1,2])
if data[0,1] == "1":
    data = np.loadtxt("DMX_vals_errs.txt", dtype=str, usecols=[1,3])

epochs = np.loadtxt("DMX_EPs.txt", dtype=float, usecols=1)
#Rs = np.loadtxt("DMX_Rs.txt", dtype=float, usecols=1)
#N = int(len(Rs)/2)
#starts = Rs[:N]
#ends   = Rs[N:]
#epochs = [(s+e)/2 for s,e in zip(starts, ends)]

vals = data[:,0]
errs = data[:,1]

vals = [float(val.replace("D","E")) for val in vals]
errs = [float(err.replace("D","E")) for err in errs]

TASC = 56195.8764880 #MJD
PB   = 0.36089348743 #days

def get_orb_phase(MJD, TASC, PB):
    diff = MJD - TASC
    if diff >= 0:
        phase = diff/PB
    else:
        phase = 1. + diff/PB
    return phase

orb_phases = [get_orb_phase(epoch, TASC, PB) for epoch in epochs]

#data = sorted(zip(epochs,vals,errs), key=lambda pair: pair[0])
data = sorted(zip(orb_phases,vals,errs), key=lambda pair: pair[0])

#vals = [val for epoch,val,err in data]
#errs = [err for epoch,val,err in data]
#epochs = [epoch for epoch,val,err in data]
vals = [val for orb_phase,val,err in data]
errs = [err for orb_phase,val,err in data]
orb_phases = [orb_phase for orb_phase,vall,err in data]

'''def f_lin(x, m, B):
    return m*x+B

def f_sin(x, A, w):
    return A*np.sin((w**2)*x)

def f_linsin(x, m, B, A, w):
    return m*x + B + A*np.sin(2*np.pi*w*(x+77))

m0 = 0.0
B0 = 0.0
A0 = 0.003
w0 = 1.0/365.0

#popt_lin, pcov_lin = curve_fit(f_lin, epochs, vals)
#popt_sin, pcov_sin = curve_fit(f_sin, epochs, vals)
popt_linsin, pcov_linsin = curve_fit(f_linsin, epochs, vals, p0=[m0,B0,A0,w0], sigma=errs, absolute_sigma=True, bounds=([-5.,-100.,0.,0.],[5.,100.,0.01,1/14.]))

slope = popt_linsin[0]
intercept = popt_linsin[1]
sin_ampl = popt_linsin[2]
sin_freq = popt_linsin[3]

xs = np.arange(min(epochs), max(epochs), 0.1)

#fit_curve = f_linsin(xs, slope, intercept, sin_ampl, sin_freq)
fit_curve = f_linsin(xs, m0, B0, A0, w0)

period = 1 / sin_freq

print(f"Sine-linear fit has slope {slope} pc/cc/day, intercept {intercept} pc/cc, sinusoidal amplitude {sin_ampl} pc/cc, and frequency {sin_freq}/day (for a period of {period} days)")

print(pcov_linsin)'''

#plt.errorbar(epochs, vals, yerr=errs, color="black", marker="x", linestyle="None", capsize=5)
#plt.plot(xs, fit_curve, c="r")
plt.errorbar(orb_phases, vals, yerr=errs, color="black", marker="x", linestyle="None", capsize=5)
#plt.xlabel("MJD")
plt.xlabel("Orbital Phase")
plt.ylabel("DMX")
plt.title("PSR J1816+4510 820 MHz, 4 TOAs per 3-min subint, 6-min DMX intervals")
plt.show()
