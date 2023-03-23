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
            self.label = 'LOFAR '+str(int(freq))+' MHz'
        else:
            self.label = 'GBT '+str(int(freq))+' MHz'

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
      'legend.fontsize'   : 11*fig_width/8.5,
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
        res, err, frq, fname = np.loadtxt(rf,dtype=str,unpack=True,usecols=[1,2,5,6])
        res = np.array([float(x)*1e3 for x in res]) #convert ms to us
        err = np.array([float(x)*1e3 for x in err]) #convert ms to us
        frq = np.array([float(x) for x in frq])
        Sband = fname == 'guppi_56452_J0636+51_0'
    else:
        res, err, frq = np.loadtxt(rf,dtype='float',unpack=True,usecols=[1,2,5])
        res,err = res*1e3, err*1e3		# Units = microseconds	
        Sband = np.full(len(res),False)
        
    obs = [freq_bw(150.0,50.0,is_Sband=False,color=colors['pink'],zorder=2), \
           freq_bw(350.0,100.0,is_Sband=False,color=colors['red'],zorder=3), \
           freq_bw(820.0,200.0,is_Sband=False,color=colors['blue'],zorder=2), \
           freq_bw(1500.0,800.0,is_Sband=False,color=colors['purple'],zorder=4), \
           freq_bw(2000.0,800.0,is_Sband=True,color=colors['yellow'],zorder=5)]

    check_tot = 0
    
    # Get orbital phases
    
    psr = pretty_psr_name.replace('$','').split()[1]
    par_fname = f'data/{psr}_tdb.par'
    tim_fname = f'data/{psr}_fiore+23.tim'
    
    if "all" in rf:
        leg_ax = ax
        tim_fname = 'data/J1816+4510_all.tim'        
    
    mo = model.get_model(par_fname)
    to = toa.get_TOAs(tim_fname,model=mo)
    mjds = [t.value for t in to.get_mjds()]
    orb_phase_unsorted = mo.orbital_phase(to, anom="mean", radians=False)
    orb_phase = np.array([op for _,op in sorted(zip(mjds,orb_phase_unsorted))])
    
    for o in obs:
        inds = o.in_bw_inds(frq)*o.Sband_check_inds(Sband)

        ax.errorbar(orb_phase[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color,ecolor=o.color, \
                    label=o.label,ms=1.1,capsize=0.75,elinewidth=0.25,zorder=o.zorder)

        check_tot += np.sum(inds)

    y_max = 1.1*max([np.abs(r)+np.abs(e) for r,e in zip(res,err)])
    
    if "all" in rf:
        y_min = 1.1*min([r-np.abs(e) for r,e in zip(res,err)])
    else:
        y_min = -y_max
    
    ax.set_ylim([y_min,y_max])
    
    # Plot a vertical line at orbital phase = 0.25
    
    ax.plot([0.25,0.25],[y_min,y_max],ls='--',color='gray',alpha=0.5,zorder=10)
    
    ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
#     if "all" in rf:
#         ax.text(0.93,0.75,"(All Residuals)",color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
    print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

fig.subplots_adjust(hspace=0)
fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)

# Add legend to the appropriate plot
leg = leg_ax.legend(numpoints=1,ncol=2,edgecolor='black',loc=6,bbox_to_anchor=(0.5,0.5),framealpha=1.0).set_zorder(20)

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