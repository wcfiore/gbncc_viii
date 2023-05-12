import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import glob
import astropy.units as u

def get_period_and_dm_from_file(infile):
    for l in infile.readlines():
        if l.startswith("F0"):
            frequency = float(l.split()[1])*u.Hz
            period = 1./frequency
            period_sec = period.to(u.s).value
            period_ms = period.to(u.ms).value
            if period < 10.*u.ms:
                period_label = f"{round(period_ms,1)} ms"
            elif period > 999.99*u.ms:
                period_label = f"{round(period_sec,2)} s"
            else:
                period_label = f"{int(round(period_ms,0))} ms"
        elif l.startswith("DM "):
            DM = float(l.split()[1])
            if DM <= 9.9:
                DM_rounded = round(DM,1)
            else:
                DM_rounded = int(round(DM,0))
            DM_label = f"{DM_rounded} pc/cm$^3$"
    return period_label, DM_label

class prof_formatting:
    def __init__(self, name, shift, dshift, offset, rotate, rotate_149, freqxshift, freqyshift, dig, freqs, dylabel, fontsize):
        self.name = name
        self.shift = shift
        self.dshift = dshift
        self.offset = offset
        self.rotate = rotate
        self.rotate_149 = rotate_149
        self.freqxshift = freqxshift
        self.freqyshift = freqyshift
        self.dig = dig
        self.freqs = freqs
        self.dylabel = dylabel
        self.fontsize = fontsize

def prof_function(file, nbin, max_prof=0.4):

    x,y,z,profile = np.loadtxt(file,dtype='float',unpack=True,skiprows=1)
    
    if len(profile)==nbin:
        pass
    else:
        scrunch_factor = len(profile)//nbin
        profile = (profile[::scrunch_factor] + profile[1::scrunch_factor]) / 2

    # Roughly identify off-pulse region
    init_median = np.median(profile)
    off_bins = np.where(profile < init_median)
    off_mean = np.mean(profile[off_bins])

    # Scale profile so max = 0.4 by default and off-pulse noise is centered on 0.
    profile /= np.max(profile)/max_prof

    return profile

def plot_profile(profile_fname, max_prof, phases, formatting, shift, freq, flux, color, nbin=128):
    y = prof_function(profile_fname, nbin=nbin, max_prof=max_prof)
    dig = formatting.dig[freq]
    dylabel  = formatting.dylabel[freq]
    
    if flux=='--':
        flux_label = False
    elif float(flux) >= 10.:
        flux_label = f"{int(float(flux))} mJy" 
    else:
        flux_label = f"{float(flux):.{dig}f} mJy"
    
    nfreq = len(formatting.freqs)
    
    if nfreq <= 2:
        ax1.plot(phases,np.roll(y,formatting.rotate)+shift,c=color)
        ax1.plot(phases,(5./4.)*np.roll(y,formatting.rotate)+shift+0.25*max_prof,c=color,alpha=0.0)
        ax1.text(formatting.freqxshift,shift+formatting.freqyshift+dylabel,flux_label, \
                 horizontalalignment='left',verticalalignment='bottom',fontsize=formatting.fontsize,color=color)
    else:
        if freq == "149 MHz":
            dylabel -= 0.04
            rotate = formatting.rotate + formatting.rotate_149
        else:
            rotate = formatting.rotate
        ax2.plot(phases,np.roll(y,rotate)+shift,c=color)
        ax2.plot(phases,np.roll(y,rotate)+shift+0.25*max_prof,c=color,alpha=0.0)
        #ax2.text(formatting.freqxshift,shift+formatting.freqyshift+dylabel+0.08,tel,horizontalalignment='left', \
        #         verticalalignment='bottom',fontsize=formatting.fontsize,color=color)
        ax2.text(formatting.freqxshift,shift+formatting.freqyshift+dylabel+0.04,freq,horizontalalignment='left', \
                 verticalalignment='bottom',fontsize=formatting.fontsize,color=color)
        if freq != "149 MHz":
            ax2.text(formatting.freqxshift,shift+formatting.freqyshift+dylabel,flux_label, \
                 horizontalalignment='left',verticalalignment='bottom',fontsize=formatting.fontsize,color=color)

# Set some plotting parameters
plt.rc("font",**{"family":"serif","serif":["Computer Modern Roman"]})
fig_width = 7.0
fig_height = fig_width
fig_size = [fig_width,fig_height]
params = {"backend": "pdf",
          "font.size"      : 14*fig_width/8.5,
          "axes.labelsize" : 14*fig_width/8.5,
          "legend.fontsize": 12*fig_width/8.5,
          "xtick.labelsize": 14*fig_width/8.5,
          "ytick.labelsize": 14*fig_width/8.5,
          "lines.linewidth": 0.8,
          "text.usetex"    : True,
          "figure.figsize" : fig_size,
          "axes.unicode_minus": True}
mpl.rcParams.update(params)
pallete = {'blue':    '#377eb8', 
          'orange':  '#ff7f00',
          'green':   '#4daf4a',
          'pink':    '#f781bf',
          'brown':   '#a65628',
          'purple':  '#984ea3',
          'gray':    '#999999',
          'red':     '#e41a1c',
          'yellow':  '#dede00'}

fig1 = plt.figure()
gs1 = gridspec.GridSpec(16, 16, figure=fig1)

# Get info from data files
DATA_PATH = "data/" # using (nearly identical) old profiles since new ones are not aligned
all_files = glob.glob(f"{DATA_PATH}*.profile")
flux_info = np.loadtxt(f"{DATA_PATH}flux.info",dtype="str")

flux57   = [ss[1] for ss in flux_info]
flux149  = [ss[2] for ss in flux_info]
flux350  = [ss[3] for ss in flux_info]
flux430  = [ss[4] for ss in flux_info]
flux820  = [ss[5] for ss in flux_info]
flux1380 = [ss[6] for ss in flux_info]
flux1500 = [ss[7] for ss in flux_info]
flux2000 = [ss[8] for ss in flux_info]
names    = [ss[0] for ss in flux_info]
        
periods = []
dms = []

Formatting_Dict = {}

all_freqs  = ["57 MHz", "149 MHz", "350 MHz", "430 MHz", "820 MHz", "1.4 GHz", \
              "1.5 GHz", "2 GHz"]

color_dict = {}
colors = [pallete['gray'], pallete['orange'], pallete['red'], pallete['pink'], pallete['blue'], pallete['purple'], pallete['purple'], pallete['yellow']]
for freq, color in zip(all_freqs, colors):
    color_dict[freq] = color

for name,s57,s149,s350,s430,s820,s1380,s1500,s2000 in zip(names,flux57,flux149,flux350,flux430, \
                                                          flux820,flux1380,flux1500,flux2000):
    
    # Get Period and DM labels
    par_path = f"{DATA_PATH}{name}_fiore+23.par"
    with open(par_path, 'r') as infile:
        period, dm = get_period_and_dm_from_file(infile)
        periods.append(period)
        dms.append(dm)
        
    # Which frequencies do we have?
    all_fluxes = [s57, s149, s350, s430, s820, s1380, s1500, s2000]
    freqs = []
    for freq,flux in zip(all_freqs, all_fluxes):
        if flux != '--':
            freqs.append(freq)
        if freq == "149 MHz" and name in ["J0214+5222","J1816+4510"]:
            freqs.append(freq)
        
    # Profile formatting
    shift  = 0.0
    dshift = 0.5
    rotate = 0.0
    rotate_149 = 0
    offset = 0.0
    showfreq = 0
    freqxshift = 0.0
    freqyshift = 0.1
    dylabel_57 = 0.0
    dylabel_149 = 0.0
    dylabel_350 = 0.0
    dylabel_430 = 0.0
    dylabel_820 = 0.0
    dylabel_1380 = 0.0
    dylabel_1500 = 0.0
    dylabel_2000 = 0.0
    dig_57 = 2
    dig_149 = 2
    dig_350 = 2
    dig_430 = 2
    dig_820 = 2
    dig_1380 = 2
    dig_1500 = 2
    dig_2000 = 2
    if name=="J0415+6111":
        dylabel_350 = 0.1
    elif name=="J0636+5128":
        rotate = 10
        freqyshift -= 0.02
    elif name=="J0957-0619":
        dylabel_350 += 0.1
    elif name=="J1239+3239":
        dylabel_820 += 0.05
    elif name=="J1327+3423":
        dig_350 = 1
        freqyshift -= 0.02
    elif name=="J1434+7257":
        rotate = 64
        dylabel_1500 = 0.05
    elif name=="J1816+4510":
        rotate = 32
        rotate_149 = 38
        freqyshift += 0.005
        dylabel_149 = 0.055
        dylabel_350 = 0.033
    elif name=="J2115+6702":
        dylabel_350 = 0.1
    elif name=="J2145+2158":
        dylabel_350 = 0.07
    elif name=="J2210+5712":
        dylabel_350 = 0.2
    elif name=="J2326+6243":
        dylabel_350 = 0.09
    elif name=="J2354-2250":
        freqyshift -= 0.05
        
    dig_dict = {} # number of decimal places to list on flux densities at each frequency
    dylabel_dict = {}
    dig_list = [dig_57, dig_149, dig_350, dig_430, dig_820, dig_1380, dig_1500, dig_2000]
    dylabels = [dylabel_57, dylabel_149, dylabel_350, dylabel_430, dylabel_820, dylabel_1380, dylabel_1500, dylabel_2000]
    
    for freq,dig,dylabel in zip(all_freqs,dig_list,dylabels):
        if freq in freqs:
            dig_dict[freq] = dig
            dylabel_dict[freq] = dylabel
    
    if len(freqs) > 2:
        fontsize = 7
        dshift = 0.2
        freqyshift -= 0.06
        offset -= 0.075
#         horizontalalignment='left'
#         verticalalignment='bottom' 
    else:
        fontsize = 6
#         horizontalalignment='left'
#         verticalalignment='bottom'
        
    Formatting_Dict[name] = prof_formatting(name=name,shift=shift,dshift=dshift,offset=offset,rotate=rotate, \
                                            rotate_149=rotate_149,freqxshift=freqxshift,freqyshift=freqyshift,dig=dig_dict, \
                                            freqs=freqs,dylabel=dylabel_dict,fontsize=fontsize)

nbin = 128
midbin = nbin/2
phases = np.arange(nbin)/float(nbin)

count = 0

# Plot profiles for pulsars with only 350 and 820 MHz

for name,period,dm,s350,s820 in zip(names,periods,dms,flux350,flux820):

    formatting = Formatting_Dict[name]
    freqs = formatting.freqs
    if len(freqs) > 2:
        continue

    profile_fname_350 = f"{DATA_PATH}{name}_350MHz_fiore+23.profile"
    profile_fname_820 = f"{DATA_PATH}{name}_820MHz_fiore+23.profile"

    print(name)
    row = int(np.floor(count/4))
    col = count % 4

    r_lo, r_hi, c_lo, c_hi = (4*row,4*(row+1), 4*col,4*(col+1))
    print(row, col)
    
    #ax1 = plt.subplot(gs1[r_lo:r_hi, c_lo:c_hi])
    ax1 = fig1.add_subplot(gs1[r_lo:r_hi, c_lo:c_hi])
    
    max_prof = 0.5
    shift = formatting.shift+formatting.offset
    
    if "350 MHz" in freqs:
        plot_profile(profile_fname=profile_fname_350, max_prof=max_prof, phases=phases, formatting=formatting, shift=shift, \
                         freq="350 MHz", flux=s350, color=color_dict["350 MHz"], nbin=nbin)
    else:
        y = prof_function(profile_fname_820, max_prof=max_prof, nbin=nbin)
        ax1.plot(phases,np.roll(y,formatting.rotate)+shift, alpha=0.0)
        
    shift += formatting.dshift
    
    if "820 MHz" in freqs:
        plot_profile(profile_fname=profile_fname_820, max_prof=max_prof, phases=phases, formatting=formatting, shift=shift, \
                     freq="820 MHz", flux=s820, color=color_dict["820 MHz"], nbin=nbin)
        
    ax1.set_ylim(-0.1,1.25)
    count += 1
    
    header_y = 1.02
    ax1.text(0.2,header_y+0.08,f"PSR {name.replace('-','$-$')}", \
             horizontalalignment='left',verticalalignment='bottom',fontsize=9)
    ax1.text(1.,header_y-0.07,dm,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
    ax1.text(0.95,header_y-0.22,period,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
    ax1.set_xticks([])
    ax1.set_yticks([])
        

# Plot profiles for pulsars with more frequencies

fig_width = 7.0
fig_height = fig_width*0.75
fig_size = [fig_width,fig_height]
params = {"backend": "pdf",
          "font.size"      : 14*fig_width/8.5,
          "axes.labelsize" : 14*fig_width/8.5,
          "legend.fontsize": 12*fig_width/8.5,
          "xtick.labelsize": 14*fig_width/8.5,
          "ytick.labelsize": 14*fig_width/8.5,
          "lines.linewidth": 0.8,
          "text.usetex"    : True,
          "figure.figsize" : fig_size,
          "axes.unicode_minus": True}
mpl.rcParams.update(params)

fig2 = plt.figure()
gs2 = gridspec.GridSpec(4, 20, figure=fig2)

count = 0

for name,period,dm,s57,s149,s350,s430,s820,s1380,s1500,s2000 in \
zip(names,periods,dms,flux57,flux149,flux350,flux430,flux820,flux1380,flux1500,flux2000):
    
    formatting = Formatting_Dict[name]
    freqs = formatting.freqs
    if len(freqs) <= 2:
        continue
    
    profile_fname_57   = f"{DATA_PATH}{name}_57.15MHz_fiore+23.profile"
    profile_fname_149  = f"{DATA_PATH}{name}_149MHz_fiore+23.profile"
    profile_fname_350  = f"{DATA_PATH}{name}_350MHz_fiore+23.profile"
    profile_fname_430  = f"{DATA_PATH}{name}_430MHz_fiore+23.profile"
    profile_fname_820  = f"{DATA_PATH}{name}_820MHz_fiore+23.profile"
    profile_fname_1380 = f"{DATA_PATH}{name}_1380MHz_fiore+23.profile"
    profile_fname_1500 = f"{DATA_PATH}{name}_1500MHz_fiore+23.profile"
    profile_fname_2000 = f"{DATA_PATH}{name}_2000MHz_fiore+23.profile"
    
    profile_fnames = [profile_fname_57,profile_fname_149,profile_fname_350,profile_fname_430,profile_fname_820, \
                      profile_fname_1380,profile_fname_1500,profile_fname_2000]
    
    fluxes = [s57,s149,s350,s430,s820,s1380,s1500,s2000]

    print(name)
    row = int(np.floor(count/5))
    col = count % 5

    r_lo, r_hi, c_lo, c_hi = (4*row,4*(row+1), 4*col,4*(col+1))
    print(row, col)

    ax2 = fig2.add_subplot(gs2[r_lo:r_hi, c_lo:c_hi])
    
    max_prof = 0.2
    
    shift = formatting.shift+formatting.offset
    
    for freq,profile_fname,flux,color in zip(all_freqs, profile_fnames, fluxes, colors):
        if freq in freqs:
            plot_profile(profile_fname=profile_fname, max_prof=max_prof, phases=phases, formatting=formatting, shift=shift, \
                     freq=freq, flux=flux, color=color_dict[freq], nbin=nbin)
        elif freq == all_freqs[-1]:
            break
        else:
            y = prof_function(profile_fname_820, max_prof=max_prof, nbin=nbin)
            ax2.plot(phases,np.roll(y,formatting.rotate)+shift, alpha=0.0)
        if freq == "57 MHz" or freq == "1.4 GHz":
            continue
        else:
            shift += formatting.dshift
            
    ax2.set_ylim([-0.1,1.25])
    count += 1

    header_y = 1.0
    ax2.text(0.015,header_y+0.2,f"PSR {name.replace('-','$-$')}", \
             horizontalalignment='left',verticalalignment='bottom',fontsize=9)
    ax2.text(1.,header_y+0.16,dm,horizontalalignment='right',verticalalignment='bottom',fontsize=7)
    ax2.text(0.95,header_y+0.12,period,horizontalalignment='right',verticalalignment='bottom',fontsize=7)
    ax2.set_xticks([])
    ax2.set_yticks([])
        
#plt.show()
fig1.savefig('profiles_1.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
fig2.savefig('profiles_2.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)