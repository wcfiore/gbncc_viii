import numpy as np
import pint.models as model
import pint.toa as toa
import pint.derived_quantities as dq
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib as mpl
import glob
    
class freq_bw:

    def __init__(self,freq,bw,is_Sband=False,color='black',zorder=1):
        self.flo      = freq - bw/2.0
        self.fhi      = freq + bw/2.0
        self.is_Sband = is_Sband
        self.color    = color
        self.zorder   = zorder
        if freq < 200.0:
            self.label = str(int(freq))+' MHz (LOFAR)'
        else:
            self.label = str(int(freq))+' MHz (GBT)'

    def in_bw_inds(self,freq_arr):
        lo_test = freq_arr > self.flo
        hi_test = freq_arr < self.fhi
        return lo_test*hi_test
    
    def Sband_check_inds(self,Sband_arr):
        if self.is_Sband:
            Sband_test = np.full(len(Sband_arr),self.is_Sband)*Sband_arr
        else:
            Sband_test = np.array([bool(1 - x) for x in Sband_arr])
        return Sband_test
    
#res_files = ['data/J1816+4510_allresiduals.dat','data/J1816+4510_residuals.dat']
res_files = ['data/J0032+6946_residuals.dat','data/J0214+5222_residuals.dat','data/J0636+5128_residuals.dat', 'data/J1239+3239_residuals.dat','data/J1816+4510_residuals.dat','data/J1816+4510_allresiduals.dat']
psr_names = [f"PSR {fname.split('/')[1].split('_')[0].replace('-','$-$')}" for fname in res_files]

nrows = len(res_files)

# Set some plotting parameters
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
fig_width = 6.0
fig_height = 6.0*nrows/6.0
fig_size = [fig_width,fig_height]
params = {'backend': 'pdf',
      'font.size'         : 14*fig_width/8.5,
      'axes.labelsize'    : 14*fig_width/8.5,
      'legend.fontsize'   : 10.5*fig_width/8.5,
      'xtick.labelsize'   : 14*fig_width/8.5,
      'ytick.labelsize'   : 14*fig_width/8.5,
      'text.usetex'       : True,
      'figure.figsize'    : fig_size,
      'axes.unicode_minus': True}
mpl.rcParams.update(params)
colors = {'blue':    '#377eb8', 
          'orange':  '#ff7f00',
          'green':   '#4daf4a',
          'pink':    '#f781bf',
          'brown':   '#a65628',
          'purple':  '#984ea3',
          'gray':    '#999999',
          'red':     '#e41a1c',
          'yellow':  '#dede00'}

fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)

for ii, (rf,nn,ax) in enumerate(zip(res_files,psr_names,axes.flat)):
    pretty_psr_name = nn 
    if "all" in rf:
        pretty_psr_name = f"{pretty_psr_name} (all residuals)"
    print(pretty_psr_name)

    if nn=='PSR J0636+5128':
            res, err, day, frq, fname = np.loadtxt(rf,dtype=str,unpack=True,usecols=[1,2,3,5,6])
            res = np.array([float(x)*1e3 for x in res]) #convert ms to us
            err = np.array([float(x)*1e3 for x in err]) #convert ms to us
            day = np.array([float(x) for x in day])
            frq = np.array([float(x) for x in frq])
            Sband = fname == 'guppi_56452_J0636+51_0'
    elif nn in ['PSR J0032+6946','PSR J0214+5222']:
        res, err, day, frq = np.loadtxt(rf,dtype='float',unpack=True) # res, err already in us
        Sband = np.full(len(res),False)
    else:
        res, err, day, frq = np.loadtxt(rf,dtype='float',unpack=True,usecols=[1,2,3,5])
        res,err = res*1e3, err*1e3  #convert ms to us
        Sband = np.full(len(res),False)
        
    obs = [freq_bw(149.0,50.0,is_Sband=False,color=colors['pink'],zorder=2), \
           freq_bw(350.0,100.0,is_Sband=False,color=colors['red'],zorder=3), \
           freq_bw(820.0,200.0,is_Sband=False,color=colors['blue'],zorder=2), \
           freq_bw(1500.0,800.0,is_Sband=False,color=colors['purple'],zorder=4), \
           freq_bw(2000.0,800.0,is_Sband=True,color=colors['yellow'],zorder=5)]

    check_tot = 0
    
    psr = pretty_psr_name.replace('$','').split()[1]
    par_fname = f'data/{psr}_fiore+23.par'
    tim_fname = f'data/{psr}_fiore+23.tim'
    
    y_max = 1.1*max([np.abs(r)+np.abs(e) for r,e in zip(res,err)])
    
    if "all" in rf:
        y_min = 1.1*min([r-np.abs(e) for r,e in zip(res,err)])
        leg_ax = ax
        tim_fname = 'data/J1816+4510_all.tim'
        # put fake points so the legend has caps on all the errorbars
        ax.errorbar(0.0,res[0],yerr=err[0],fmt='o',mfc=colors['orange'],mec=colors['orange'],ecolor=colors['orange'], \
                        label="149 MHz (LOFAR)",ms=1.5,capsize=1.,elinewidth=0.5)
        ax.errorbar(0.0,res[0],yerr=err[0],fmt='o',mfc=colors['yellow'],mec=colors['yellow'],ecolor=colors['yellow'], \
                        label="2000 MHz (GBT)",ms=1.5,capsize=1.,elinewidth=0.5)
    else:
        y_min = -y_max
        
    ax.set_ylim([y_min,y_max])
    
    # Get orbital phases
    
    mo = model.get_model(par_fname)
    OM = mo["OM"].value
    to = toa.get_TOAs(tim_fname,model=mo)
    mjds = [t.value for t in to.get_mjds()]
    orb_phase_unsorted = mo.orbital_phase(to, anom="mean", radians=False) # I think this is not actually the mean anomaly M,
                                                                 # but the Phi used in the ELL1 model,
                                                                 # Phi = M + omega
    orb_phase = np.array([op for _,op in sorted(zip(mjds,orb_phase_unsorted))])
    
    for o in obs:
        inds = o.in_bw_inds(frq)*o.Sband_check_inds(Sband)
        ax.errorbar(orb_phase[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color,ecolor=o.color, \
                    label=o.label,ms=1.1,capsize=0.75,elinewidth=0.25,zorder=o.zorder)

        check_tot += np.sum(inds)
        
    # Plot a vertical line at orbital phase = 0.25
    
    ax.plot([0.25,0.25],[y_min,y_max],ls='--',color='gray',alpha=0.5,zorder=10)
    
    ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
#     if "all" in rf:
#         ax.text(0.93,0.75,"(All Residuals)",color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
    print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

fig.subplots_adjust(hspace=0)
fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)

# once again assuring the legend errorbars have caps here
handles, labels = leg_ax.get_legend_handles_labels()
by_label = dict(zip(reversed(labels), reversed(handles)))

# Add legend to the appropriate plot
leg = leg_ax.legend(reversed(by_label.values()),reversed(by_label.keys()),numpoints=1,ncol=2,edgecolor='black',loc=6, \
                    bbox_to_anchor=(0.485,0.5),framealpha=1.0).set_zorder(20)

# Calculate *global* xlims
x_lims = [0.0,1.0]
ax.set_xlabel('Orbital Phase')
ax.set_xlim(x_lims)

# Plot center lines using global xlims
for ax in axes.flat:
    ax.plot(x_lims,[0,0],ls='--',color='black',alpha=0.5,zorder=10)
    ax.ticklabel_format(useOffset=False,style='plain')

#######################

# Add ylabel, adjust layout and save figure
fig.text(-0.02, 0.5, r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')

fig.savefig('res_orb.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)