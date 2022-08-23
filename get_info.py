import numpy as np
import matplotlib.pyplot as plt

class orb_data():
  def __init__(self,fname,tolerance=0.001):
    data = np.loadtxt(fname,dtype='float')
    self.ep = data[:,0]
    self.p0 = data[:,1]
    self.t0 = self.ep[0]
    self.tol = tolerance 

  def R(self,pb_guess):
    phase = ((self.ep-self.t0)/pb_guess)%1.0 
    inds = np.argsort(phase)
    sort_p0 = self.p0[inds]
    R = 0
    for i in range(len(sort_p0)-1):
      R += (sort_p0[i]-sort_p0[i+1])**2
    R += (sort_p0[0]-sort_p0[-1])**2   # Account for one more term (wrap)
    return R

  def pb_inc(self,pb,T):
    return self.tol*pb**2/(2.*np.pi*T)

  def get_pb(self,pbrange=None):
    epsort = np.sort(self.ep)
    T      = np.max(self.ep)-np.min(self.ep)    # Total data span

    if not pbrange:
      PBlo   = (epsort[1]-epsort[0])/2.0
      PBhi   = T*2.0
    else:
      PBlo,PBhi = pbrange

    print(f"Observations span {T:.2f} days.")

    pb_arr = [PBlo]
    R_arr  = [self.R(PBlo)]
    pb = PBlo+self.pb_inc(PBlo,T)
    while pb < PBhi:
      pb_arr.append(pb)
      R_arr.append(self.R(pb))
      pb += self.pb_inc(pb,T)
    pbind = np.argmin(R_arr)
    print(f"Average spin period: {np.mean(self.p0):.4f} ms.")
    print(f"Minimum roughness of {min(R_arr)} with PB = {pb_arr[pbind]}.")
    plt.plot(pb_arr,R_arr,color='black')
    plt.xlim([PBlo,PBhi])
    plt.xlabel(r'Trial Orbital Period (days)')
    plt.ylabel(r'Roughness (ms$^2$)')
    plt.show()
    #plt.savefig('rough.pdf',format='pdf') 
    return pb_arr[pbind]

  def est_a1(self,pb_guess,ecc_guess):
    p0_extent = (np.max(self.p0)-np.min(self.p0))/1000.
    p0_guess = np.mean(self.p0)/1000.
    c = 3.0e8
    pb_s = pb_guess * 24.0 * 3600.0
    orb_freq = 2*np.pi/pb_s
    print((p0_extent*np.sqrt(1-ecc_guess**2))/(2*orb_freq*p0_guess))

the_data = orb_data('data/1017-mjd-p0-err.txt')
pb_best = the_data.get_pb(pbrange=[5,15.])
the_data.est_a1(pb_best,0.2)
