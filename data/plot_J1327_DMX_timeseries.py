import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("DMX_vals_errs.txt", dtype=str, usecols=[1,2])
if data[0,1] == "1":
    data = np.loadtxt("DMX_vals_errs.txt", dtype=str, usecols=[1,3])

epochs = np.loadtxt("DMX_EPs.txt", dtype=float, usecols=1)

vals = data[:,0]
errs = data[:,1]

vals = [float(val.replace('D','E')) for val in vals]
errs = [float(err.replace('D','E')) for err in errs]

data = sorted(zip(epochs,vals,errs), key=lambda pair: pair[0])

vals = [val for epoch,val,err in data]
errs = [err for epoch,val,err in data]
epochs = [epoch for epoch,vall,err in data]

plt.errorbar(epochs, vals, yerr=errs, color="black", marker="x", linestyle='None', capsize=5)
plt.xlabel("MJD")
plt.ylabel("DMX")
plt.xlim(58250,max(epochs)+50)
plt.title("PSR J1327+3423 6.5-day Intervals")
plt.show()
