import numpy as np
import pint.models as model
import pint.derived_quantities as dq
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib as mpl
import glob

# Set some plotting parameters
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
fig_width = 6.0
fig_height = 8.0  #fig_width*0.6
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

res_files1 = ['J1327+3423_residuals.dat']
names1 = ['PSR J1327+3423']
res_files2 = ['J0032+6946_residuals.dat','J0214+5222_residuals.dat','J0636+5128_residuals.dat', 'J1239+3239_residuals.dat','J1434+7257_residuals.dat','J1816+4510_residuals.dat']
names2 = ['PSR J0032+6946','PSR J0214+5222','PSR J0636+5128','PSR J1239+3239','PSR J1434+7257','PSR J1816+4510']
res_files3 = ['J0141+6303_residuals.dat','J0415+6111_residuals.dat','J0957-0619_residuals.dat','J1505-2524_residuals.dat','J1530-2114_residuals.dat','J1913+3732_residuals.dat']
names3 = ['PSR J0141+6303','PSR J0415+6111','PSR J0957$-$0619','PSR J1505$-$2524','PSR J1530$-$2114','PSR J1913+3732']
res_files4 = ['J1929+6630_residuals.dat','J1930+6205_residuals.dat','J2104+2830_residuals.dat','J2115+6702_residuals.dat','J2145+2158_residuals.dat','J2210+5712_residuals.dat']
names4 = ['PSR J1929+6630','PSR J1930+6205','PSR J2104+2830','PSR J2115+6702','PSR J2145+2158','PSR J2210+5712']
res_files5 = ['J2326+6243_residuals.dat','J2354-2250_residuals.dat']
names5 = ['PSR J2326+6243','PSR J2354$-$2250']

for jj in range(5):

    plotnum = str(jj+1)
    if jj == 0:
        res_files = ['data/'+rf for rf in res_files1]
        names = names1
    elif jj == 1:
        res_files = ['data/'+rf for rf in res_files2]
        names = names2
    elif jj == 2:
        res_files = ['data/'+rf for rf in res_files3]
        names = names3
    elif jj == 3:
        res_files = ['data/'+rf for rf in res_files4]
        names = names4
    else:
        res_files = ['data/'+rf for rf in res_files5]
        names = names5

    nrows = len(res_files)

    # Make five plots (in batches of 5/6, except for J1327)
    if nrows==1:
        fig = plt.figure(figsize=(6,2))
        ax = fig.add_subplot(1, 1, 1)
        obs = [freq_bw_rcvr(57.15,80.0,rcvr='LWA1',color='black',zorder=2), \
               freq_bw_rcvr(327.0,90.0,rcvr='puppi_327',color='purple',zorder=5), \
               freq_bw_rcvr(350.0,100.0,rcvr='guppi',color='r',zorder=4), \
               freq_bw_rcvr(430.0,100.0,rcvr='puppi_430',color='pink',zorder=3), \
               freq_bw_rcvr(820.0,200.0,rcvr='guppi',color='b',zorder=3), \
               freq_bw_rcvr(1380.0,800.0,rcvr='puppi_1380',color='cyan',zorder=2)]
        min_yr,max_yr = (9999.0,0.0)
        pretty_psr_name = names[0] 
        print(pretty_psr_name)
        rf = res_files[0]
        res, err, day, frq, tel = np.loadtxt(rf,dtype=str,unpack=True,usecols=[1,2,3,5,7]) 
        res = np.array([float(x)*1e3 for x in res]) #convert ms to us
        err = np.array([float(x)*1e3 for x in err]) #convert ms to us
        day = np.array([float(x) for x in day])
        frq = np.array([float(x) for x in frq])
        check_tot = 0
        mjd_T = Time(day,format='mjd')
        yr_T = Time(mjd_T,format='mjd').decimalyear

        if min(yr_T)<min_yr:
            min_yr = min(yr_T)
        if max(yr_T)>max_yr:
            max_yr = max(yr_T)

        # Get phase info
        psr = pretty_psr_name.replace('$','').split()[1]
        par_fname = f'data/{psr}_fiore+22.par'

        with open(par_fname, 'r') as infile:
            for l in infile.readlines():
                if l.startswith("F0"):
                    f0      = float(l.split()[1])*u.Hz
                    f0err  = float(l.split()[3])*u.Hz
                elif l.startswith("F1"):
                    f1      = float(l.split()[1])*u.Hz/u.s
                    f1err  = float(l.split()[3])*u.Hz/u.s
                else:
                    pass

        p0,p0err,p1,p1err = dq.pferrs(f0,f0err,f1,f1err)

        ax_R = ax.twinx()
#         plotted_day = np.array([])
#         plotted_frq = np.array([])
        for o in obs:
            inds = o.in_bw_inds(frq)*o.same_rcvr_inds(tel)
#             plotted_day = np.append(plotted_day,day[inds])
#             plotted_frq = np.append(plotted_frq,frq[inds])
            ax.errorbar(yr_T[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color, \
                        ecolor=o.color,label=o.label,ms=1.5,capsize=1,elinewidth=0.5,zorder=o.zorder)
            # ax_R.errorbar(day,res/p0.value/1.0e6,yerr=err,fmt='o',alpha=0.0) # res in microseconds
            check_tot += np.sum(inds)
#         for d,f in zip(day,frq):
#             if d not in plotted_day:
#                 print(f"{d} {f}")
        #rms = np.sqrt(np.mean(res**2))
        #ax.set_ylim([-3*rms,3*rms])
        #ax_R.set_ylim([-3*rms/p0.value/1.0e6,3*rms/p0.value/1.0e6])
        #ylim = 1.3*max(np.abs(res))
        ylim = 1.1*max([np.abs(r)+np.abs(e) for r,e in zip(res,err)])
        ax.set_ylim([-ylim,ylim])
        ax_R.set_ylim([-ylim/p0.to(u.us).value,ylim/p0.to(u.us).value])
        ax.minorticks_on()
        ax_R.minorticks_on()

        # 3 labeled ticks: https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
        ax_R.yaxis.set_major_locator(plt.MaxNLocator(3))
        #ax_R.ticklabel_format(useOffset=False,style='plain')

        ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
        print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

        # Calculate xlims
        x_lims = [min_yr,max_yr]
        timespan = x_lims[1]-x_lims[0]
        x_lims[0] -= 0.05 * timespan
        x_lims[1] += 0.05 * timespan
        ax.set_xlabel('Year')
        ax.set_xlim(x_lims)

        # Make lims Time objects (use for twin axis)
        minmax_yr_T = [Time(xl,format='decimalyear') for xl in x_lims]
        minmax_mjd_T = [mmyT.mjd for mmyT in minmax_yr_T] 

        print(f'Observations plotted here span years {minmax_yr_T[0].value:.2f}-{minmax_yr_T[1].value:.2f},')
        print(f'                       ...  and MJDs {minmax_mjd_T[0]:.2f}-{minmax_mjd_T[1]:.2f}.')

        # Plot center line using xlim
        ax.plot(x_lims,[0,0],ls='--',color='black',alpha=0.5,zorder=7)
        ax.ticklabel_format(useOffset=False,style='plain')
        
        # Add legend
        leg = ax.legend(numpoints=1,edgecolor='black',loc=6,bbox_to_anchor=(0.05,0.5),fontsize='x-small',framealpha=1.0).set_zorder(20)

        # Add MJD scale
        ax2 = ax.twiny()
        ax2.errorbar(day,res,yerr=err,fmt='o',alpha=0.0)
        ax2.set_xlabel('Modified Julian Date')
        ax2.set_xlim(minmax_mjd_T)
        ax2.minorticks_on()

        # Add ylabel
        fig.text(0.0, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
        fig.text(1.0, 0.5,r'Phase', ha='center', va='center', rotation='vertical')

        # Adjust layout and save figure
        fig.subplots_adjust(hspace=0)
        fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)
        fig.savefig('res'+plotnum+'.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
    else:
        fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)

        min_yr,max_yr = (9999.0,0.0)
        for ii, (rf,nn,ax) in enumerate(zip(res_files,names,axes.flat)):

            #file_psr_name = ('psrJ%s' % (rf.split('_')[0]))
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
            if jj==1:
                obs = [freq_bw(150.0,50.0,is_Sband=False,color='black',zorder=2), \
                       freq_bw(350.0,100.0,is_Sband=False,color='r',zorder=3), \
                       freq_bw(820.0,200.0,is_Sband=False,color='b',zorder=2), \
                       freq_bw(1500.0,800.0,is_Sband=False,color='cyan',zorder=4), \
                       freq_bw(2000.0,800.0,is_Sband=True,color='g',zorder=5)]
            else:
                obs = [freq_bw(350.0,100.0,color='r'),freq_bw(820.0,200.0,color='b')]

            check_tot = 0
            mjd_T = Time(day,format='mjd')
            yr_T = Time(mjd_T,format='mjd').decimalyear 

            if ii == 0:
                top_ax = ax

            if min(yr_T)<min_yr:
                min_yr = min(yr_T)
            if max(yr_T)>max_yr:
                max_yr = max(yr_T)

            # Get phase info
            psr = pretty_psr_name.replace('$','').split()[1]
            par_fname = f'data/{psr}_fiore+22.par'

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
#             plotted_day = np.array([])
#             plotted_frq = np.array([])
            for o in obs:
                inds = o.in_bw_inds(frq)*o.Sband_check_inds(Sband)
#                 plotted_day = np.append(plotted_day,day[inds])
#                 plotted_frq = np.append(plotted_frq,frq[inds])
                ax.errorbar(yr_T[inds],res[inds],yerr=err[inds],fmt='o',mfc=o.color,mec=o.color,ecolor=o.color, \
                            label=o.label,ms=1.5,capsize=1.,elinewidth=0.5,zorder=o.zorder)
                # ax_R.errorbar(day,res/p0.value/1.0e6,yerr=err,fmt='o',alpha=0.0) # res in microseconds
                check_tot += np.sum(inds)
#             for d,f in zip(day,frq):
#                 if d not in plotted_day:
#                     print(f"{d} {f}")
            #rms = np.sqrt(np.mean(res**2))
            #ax.set_ylim([-5*rms,5*rms])
            #ax_R.set_ylim([-5*rms/p0.value/1.0e6,5*rms/p0.value/1.0e6])
            #ylim = 1.3*max(np.abs(res))
            ylim = 1.1*max([np.abs(r)+np.abs(e) for r,e in zip(res,err)])
            ax.set_ylim([-ylim,ylim])
            ax_R.set_ylim([-ylim/p0.to(u.us).value,ylim/p0.to(u.us).value])
            ax.minorticks_on()
            ax_R.minorticks_on()

            # 3 labeled ticks: https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
            ax_R.yaxis.set_major_locator(plt.MaxNLocator(3))
            ax_R.ticklabel_format(useOffset=False,style='plain')
            if jj==2 or jj==3:
                ax.text(0.34,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
            else:
                ax.text(0.93,0.9,pretty_psr_name,color='black',rotation=0,size=10,va='center',ha='right',transform=ax.transAxes)
            print(str(check_tot)+'/'+str(len(res))+' residuals plotted.')

        # Add legend to the appropriate plot
        if jj==1:
            #leg = axes.flat[2].legend(numpoints=1,edgecolor='black',loc=4,framealpha=1.0).set_zorder(20)
            leg = top_ax.legend(numpoints=1,edgecolor='black',loc=6,bbox_to_anchor=(0.05,0.5),framealpha=1.0).set_zorder(20)
        elif jj==2 or jj==3:
            leg = top_ax.legend(numpoints=1,edgecolor='black',loc=3,framealpha=1.0).set_zorder(20)
        else:
            leg = top_ax.legend(numpoints=1,edgecolor='black',loc=4,framealpha=1.0).set_zorder(20)
        
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
        if jj==4:
            bot = 0.54
            labels=0.72
            fig.text(0.01, labels,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
            fig.text(1.01, labels,r'Phase', ha='center', va='center', rotation='vertical')
            fig.subplots_adjust(bottom=bot)
        else:
            fig.text(-0.007, 0.5,r'Residual ($\mu$s)', ha='center', va='center', rotation='vertical')
            fig.text(1.025, 0.5,r'Phase', ha='center', va='center', rotation='vertical')
            fig.subplots_adjust(hspace=0)
            fig.tight_layout(pad=0.0, w_pad=0.5, h_pad=0.25)
        fig.savefig('res'+plotnum+'.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)