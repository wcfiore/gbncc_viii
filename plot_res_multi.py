import numpy as np
import pint.models as model
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
    
class freq_bw_rcvr:

    def __init__(self,freq,bw,rcvr,color='black',zorder=1):
        self.flo    = freq - bw/2.0
        self.fhi    = freq + bw/2.0
        self.rcvr   = rcvr
        self.color  = color
        self.zorder = zorder
        if 'guppi' in rcvr:
            self.label = 'GBT '+str(int(freq))+' MHz'
        elif 'LWA1' in rcvr:
            self.label = 'LWA1 '+str(int(freq))+' MHz'
        elif 'puppi' in rcvr:
            self.label = 'AO '+str(int(freq))+' MHz'
        else:
            self.label = str(int(freq))+' MHz'

    def in_bw_inds(self,freq_arr):
        lo_test = freq_arr > self.flo
        hi_test = freq_arr < self.fhi
        return lo_test*hi_test
    
    def same_rcvr_inds(self,rcvr_arr):
        return self.rcvr == rcvr_arr

res_files1 = ['J0032+6946_residuals.dat','J0214+5222_residuals.dat','J0636+5128_residuals.dat', 'J1239+3239_residuals.dat','J1434+7257_residuals.dat','J1816+4510_residuals.dat']
res_files2 = ['J0141+6303_residuals.dat','J0415+6111_residuals.dat','J0957-0619_residuals.dat','J1327+3423_residuals.dat','J1505-2524_residuals.dat','J1530-2114_residuals.dat']
res_files3 = ['J1913+3732_residuals.dat','J1929+6630_residuals.dat','J1930+6205_residuals.dat','J2104+2830_residuals.dat','J2115+6702_residuals.dat']
res_files4 = ['J2145+2158_residuals.dat','J2210+5712_residuals.dat','J2326+6243_residuals.dat','J2354-2250_residuals.dat']

for jj in range(4):

    plotnum = str(jj+1)
    if jj == 0:
        res_files = ['data/'+rf for rf in res_files1]
    elif jj == 1:
        res_files = ['data/'+rf for rf in res_files2]
    elif jj == 2:
        res_files = ['data/'+rf for rf in res_files3]
    else:
        res_files = ['data/'+rf for rf in res_files4]
        
    names = [f"PSR {fname.split('/')[1].split('_')[0].replace('-','$-$')}" for fname in res_files]

    nrows = len(res_files)
    
    # Set some plotting parameters
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    fig_width = 6.0
    fig_height = 8.0*(nrows/6.)  #fig_width*0.6
    fig_size = [fig_width,fig_height]
    params = {'backend': 'pdf',
              'font.size'         : 14*fig_width/8.5,
              'axes.labelsize'    : 14*fig_width/8.5,
              'legend.fontsize'   : 12*fig_width/8.5,
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

    min_yr,max_yr = (9999.0,0.0)
    for ii, (rf,nn,ax) in enumerate(zip(res_files,names,axes.flat)):
        pretty_psr_name = nn 
        print(pretty_psr_name)

        if nn=='PSR J0636+5128':
            res, err, day, frq, fname = np.loadtxt(rf,dtype=str,unpack=True,usecols=[1,2,3,5,6])
            res = np.array([float(x)*1e3 for x in res]) #convert ms to us
            err = np.array([float(x)*1e3 for x in err]) #convert ms to us
            day = np.array([float(x) for x in day])
            frq = np.array([float(x) for x in frq])
            Sband = fname == 'guppi_56452_J0636+51_0'
        else:
            res, err, day, frq = np.loadtxt(rf,dtype='float',unpack=True,usecols=[1,2,3,5])
            res,err = res*1e3, err*1e3		# Units = microseconds	
            Sband = np.full(len(res),False)
            
        if jj==0:
            obs = [freq_bw(150.0,50.0,is_Sband=False,color=colors['orange'],zorder=2), \
                   freq_bw(350.0,100.0,is_Sband=False,color=colors['red'],zorder=3), \
                   freq_bw(820.0,200.0,is_Sband=False,color=colors['blue'],zorder=2), \
                   freq_bw(1500.0,800.0,is_Sband=False,color=colors['purple'],zorder=4), \
                   freq_bw(2000.0,800.0,is_Sband=True,color=colors['yellow'],zorder=5)]
        elif jj==1:
            obs = [freq_bw_rcvr(57.15,80.0,rcvr='LWA1',color=colors['gray'],zorder=2), \
                freq_bw_rcvr(350.0,100.0,rcvr='guppi',color=colors['red'],zorder=4), \
                freq_bw_rcvr(430.0,100.0,rcvr='puppi_430',color=colors['pink'],zorder=3), \
                freq_bw_rcvr(820.0,200.0,rcvr='guppi',color=colors['blue'],zorder=3), \
                freq_bw_rcvr(1380.0,800.0,rcvr='puppi_1380',color=colors['purple'],zorder=2)]
        else:
            obs = [freq_bw(350.0,100.0,color=colors['red']),freq_bw(820.0,200.0,color=colors['blue'])]
#         else:
#             if nn=="PSR J1327+3423":
#                 obs = [freq_bw_rcvr(57.15,80.0,rcvr='LWA1',color='black',zorder=2), \
#                 freq_bw_rcvr(350.0,100.0,rcvr='guppi',color='r',zorder=4), \
#                 freq_bw_rcvr(430.0,100.0,rcvr='puppi_430',color='pink',zorder=3), \
#                 freq_bw_rcvr(820.0,200.0,rcvr='guppi',color='b',zorder=3), \
#                 freq_bw_rcvr(1380.0,800.0,rcvr='puppi_1380',color='cyan',zorder=2)]
#             else:
#                 obs = [freq_bw(350.0,100.0,color='r'),freq_bw(820.0,200.0,color='b')]

        check_tot = 0
        mjd_T = Time(day,format='mjd')
        yr_T = Time(mjd_T,format='mjd').decimalyear 

        if ii == 0:
            leg_ax = ax
            top_ax = ax
            
#         if nn=="PSR J1327+3423":
#             leg_ax = ax

        if min(yr_T)<min_yr:
            min_yr = min(yr_T)
        if max(yr_T)>max_yr:
            max_yr = max(yr_T)

        # Get phase info
        psr = pretty_psr_name.replace('$','').split()[1]
        par_fname = f'data/{psr}_fiore+23.par'

        with open(par_fname, 'r') as infile:
            for l in infile.readlines():
                if l.startswith("F0"):
                    f0      = float(l.split()[1])*u.Hz
                    f0err   = float(l.split()[3])*u.Hz
                elif l.startswith("F1"):
                    f1      = float(l.split()[1])*u.Hz/u.s
                    f1err   = float(l.split()[3])*u.Hz/u.s
                else:
                    pass

        p0,p0err,p1,p1err = dq.pferrs(f0,f0err,f1,f1err)

        ax_R = ax.twinx()
        for o in obs:
            if jj==0:
                inds = o.in_bw_inds(frq)*o.Sband_check_inds(Sband)
            else:
                inds = o.in_bw_inds(frq)
            ax.errorbar(yr_T[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color,ecolor=o.color, \
                        label=o.label,ms=1.5,capsize=1.,elinewidth=0.5,zorder=o.zorder)
            check_tot += np.sum(inds)
            
        ylim = 1.1*max([np.abs(r)+np.abs(e) for r,e in zip(res,err)])
        ax.set_ylim([-ylim,ylim])
        ax_R.set_ylim([-ylim/p0.to(u.us).value,ylim/p0.to(u.us).value])
        ax.minorticks_on()
        ax_R.minorticks_on()

        # 3 labeled ticks: https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        ax_R.yaxis.set_major_locator(plt.MaxNLocator(3))
        ax_R.ticklabel_format(useOffset=False,style='plain')
        ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
        print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

    if jj==0 or jj==1:
        leg = leg_ax.legend(numpoints=1,edgecolor='black',loc=4,framealpha=1.0).set_zorder(20)
    else:
        leg = leg_ax.legend(numpoints=1,edgecolor='black',loc=4,bbox_to_anchor=(1.0,0.01),framealpha=1.0).set_zorder(20)

    # Calculate *global* xlims
    #x_lims = [min_yr,max_yr]
    alt = True # different time ranges for each plot
    if jj==0:
        x_lims = [2009.75,2022.75]
    elif jj==1:
        x_lims = [2012.0,2021.0]
    elif jj==2:
        x_lims = [2012.0,2017.75]
    elif jj==3:
        x_lims = [2013.0,2023.0]
#     alt = False # make each plot use the same time range
#     x_lims = [2009.75,2022.75]
    ax.set_xlabel('Year')
    ax.set_xlim(x_lims)

    # Make lims Time objects (use for twin axis)
    minmax_yr_T = [Time(xl,format='decimalyear') for xl in x_lims]
    minmax_mjd_T = [mmyT.mjd for mmyT in minmax_yr_T] 

    print(f'Observations plotted here span years {minmax_yr_T[0].value:.2f}-{minmax_yr_T[1].value:.2f},')
    print(f'                       ...  and MJDs {minmax_mjd_T[0]:.2f}-{minmax_mjd_T[1]:.2f}.')

    # Plot center lines using global xlims
    for ax in axes.flat:
        ax.plot(x_lims,[0,0],ls='--',color='black',alpha=0.5,zorder=10)
        ax.ticklabel_format(useOffset=False,style='plain')

    # Add MJD scale
    ax2 = top_ax.twiny()
    ax2.errorbar(day,res,yerr=err,fmt='o',alpha=0.0)
    ax2.set_xlabel('Modified Julian Date')
    ax2.set_xlim(minmax_mjd_T)
    ax2.minorticks_on()

    # Add common ylabels, adjust layout and save figure
    if jj==0:
        fig.text(-0.01, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
        fig.text(1.025, 0.5,r'Phase', ha='center', va='center', rotation='vertical')
        fig.subplots_adjust(hspace=0)
        fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)
    elif jj==1:
        fig.text(-0.02, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
        fig.text(1.025, 0.5,r'Phase', ha='center', va='center', rotation='vertical')
        fig.subplots_adjust(hspace=0)
        fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)
    elif jj==2:
        fig.text(0.02, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
        fig.text(1.025, 0.5,r'Phase', ha='center', va='center', rotation='vertical')
    elif jj==3:
        fig.text(0.02, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
        fig.text(1.025, 0.5,r'Phase', ha='center', va='center', rotation='vertical')
    if alt:
        fig.savefig('res'+plotnum+'_alt.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
    else:
        fig.savefig('res'+plotnum+'.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)