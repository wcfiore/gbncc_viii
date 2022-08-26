import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib as mpl
import glob

# Set some plotting parameters
plt.rc("font",**{"family":"serif","serif":["Computer Modern Roman"]})
fig_width = 6.0
fig_height = 8.0  #fig_width*0.6
fig_size = [fig_width,fig_height]
params = {"backend": "pdf",
          "font.size"      : 14*fig_width/8.5,
          "axes.labelsize" : 14*fig_width/8.5,
          "legend.fontsize": 12*fig_width/8.5,
          "xtick.labelsize": 14*fig_width/8.5,
          "ytick.labelsize": 14*fig_width/8.5,
          "text.usetex"    : True,
          "figure.figsize" : fig_size,
          "axes.unicode_minus": True}
mpl.rcParams.update(params)

class freq_bw:

    def __init__(self,freq,bw,color='black'):
        self.flo   = freq - bw/2.0
        self.fhi   = freq + bw/2.0
        self.color = color
        self.label = str(int(freq))+' MHz'

    def in_bw_inds(self,freq_arr):
        lo_test = freq_arr > self.flo
        hi_test = freq_arr < self.fhi
        return lo_test*hi_test

res_files1 = ['0405+3347_resids.dat','1018-1523_resids.dat','1045-0436_resids.dat','1122-3546_resids.dat','2022+2534_resids.dat','2039-3616_resids.dat']
names1 = ['PSR J0405+3347','PSR J1018$-$1523','PSR J1045$-$0436','PSR J1122$-$3546','PSR J2022+2534','PSR J2039$-$3616']
res_files2 = ['1221-0633_resids.dat','1317-0157_resids.dat','1742-0203_resids.dat','2017-2737_resids.dat','2018-0414_resids.dat']
names2 = ['PSR J1221$-$0633','PSR J1317$-$0157','PSR J1742$-$0203','PSR J2017$-$2737','PSR J2018$-$0414']

for jj in range(2):

    plotnum = str(jj+1)
    if jj == 0:
        res_files = ['data/'+rf for rf in res_files1]
        names = names1
    else:
        res_files = ['data/'+rf for rf in res_files2]
        names = names2

    nrows = len(res_files)

    # Make two plots (in batches of 5/6)
    fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)

    min_yr,max_yr = (9999.0,0.0)
    for ii, (rf,nn,ax) in enumerate(zip(res_files,names,axes.flat)):

        #file_psr_name = ('psrJ%s' % (rf.split('_')[0]))
        pretty_psr_name = nn 
        print(pretty_psr_name)

        day, res, err, frq = np.loadtxt(rf,dtype='float',unpack=True) 
        res,err = res*1e6, err*1e6		# Units = microseconds	
        obs = [freq_bw(820.0,200.0,color='b'),freq_bw(350.0,100.0,'r')]

        check_tot = 0
        mjd_T = Time(day,format='mjd')
        yr_T = Time(mjd_T,format='mjd').decimalyear 

        if ii == 0:
            top_ax = ax

        if min(yr_T)<min_yr:
            min_yr = min(yr_T)
        if max(yr_T)>max_yr:
            max_yr = max(yr_T)

        for o in obs:
            inds = o.in_bw_inds(frq)
            ax.errorbar(yr_T[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color,ecolor=o.color,label=o.label,ms=4,capsize=3)
            check_tot += np.sum(inds)

        rms = np.sqrt(np.mean(res**2))
        ax.set_ylim([-5*rms,5*rms])
        ax.minorticks_on()

        ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
        print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

    # Add legend to bottom plot
    leg = top_ax.legend(numpoints=1,edgecolor='black',loc=4)
    leg.get_frame().set_alpha(1.0)

    # Calculate *global* xlims
    x_lims = [min_yr,max_yr]
    timespan = x_lims[1]-x_lims[0]
    x_lims[0] -= 0.05 * timespan
    x_lims[1] += 0.05 * timespan
    ax.set_xlabel('Year')
    ax.set_xlim(x_lims)

    # Make lims Time objects (use for twin axis)
    minmax_yr_T = [Time(xl,format='decimalyear') for xl in x_lims]
    minmax_mjd_T = [mmyT.mjd for mmyT in minmax_yr_T] 

    print(f"Observations plotted here span years {minmax_yr_T[0].value:.2f}-{minmax_yr_T[1].value:.2f},")
    print(f"                       ...  and MJDs {minmax_mjd_T[0]:.2f}-{minmax_mjd_T[1]:.2f}.")

    # Plot center lines using global xlims
    for ax in axes.flat:
        ax.plot(x_lims,[0,0],ls='--',color='black',alpha=0.5)

    # Add MJD scale
    ax2 = top_ax.twiny()
    ax2.errorbar(day,res,yerr=err,fmt='o',alpha=0.0)
    ax2.set_xlabel('Modified Julian Date')
    ax2.set_xlim(minmax_mjd_T)
    ax2.minorticks_on()

    # Add common ylabel
    fig.text(0.0, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')

    # Adjust layout and save figure
    fig.subplots_adjust(hspace=0)
    fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)
    fig.savefig('res'+plotnum+'.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
