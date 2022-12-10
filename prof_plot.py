import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import glob
import astropy.units as u

def prof_function(file):

    x,y,z,profile = np.loadtxt(file,dtype='float',unpack=True,skiprows=1)

    # Roughly identify off-pulse region
    init_median = np.median(profile)
    off_bins = np.where(profile < init_median)
    off_mean = np.mean(profile[off_bins])

    # Scale profile so max = 0.5 and off-pulse noise is centered on 0.
    profile -= off_mean
    profile /= 2 * np.max(profile) 

    return profile

# Set some plotting parameters
plt.rc("font",**{"family":"serif","serif":["Computer Modern Roman"]})
fig_width = 6.0
fig_height = fig_width*0.6
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

# Get info from data files
DATA_PATH = "data/" # using (nearly identical) old profiles since new ones are not aligned
all_files = glob.glob(f"{DATA_PATH}*.profile")
flux_info = np.loadtxt(f"{DATA_PATH}flux.info",dtype="str")

for i in range(3):
    if i==0:
        flux57   = [ss[1] for ss in flux_info][:-8][:7]+[ss[1] for ss in flux_info][:-8][8:]
        flux327  = [ss[2] for ss in flux_info][:-8][:7]+[ss[2] for ss in flux_info][:-8][8:]
        flux350  = [ss[3] for ss in flux_info][:-8][:7]+[ss[3] for ss in flux_info][:-8][8:]
        flux430  = [ss[4] for ss in flux_info][:-8][:7]+[ss[4] for ss in flux_info][:-8][8:]
        flux820  = [ss[5] for ss in flux_info][:-8][:7]+[ss[5] for ss in flux_info][:-8][8:]
        flux1380 = [ss[6] for ss in flux_info][:-8][:7]+[ss[6] for ss in flux_info][:-8][8:]
        flux1500 = [ss[7] for ss in flux_info][:-8][:7]+[ss[7] for ss in flux_info][:-8][8:]
        flux2000 = [ss[8] for ss in flux_info][:-8][:7]+[ss[8] for ss in flux_info][:-8][8:]
        names    = [ss[0] for ss in flux_info][:-8][:7]+[ss[0] for ss in flux_info][:-8][8:]
        fig1 = plt.figure()
        gs1 = gridspec.GridSpec(12, 16)
    if i==1:
        flux57   = [ss[1] for ss in flux_info][13:]
        flux327  = [ss[2] for ss in flux_info][13:]
        flux350  = [ss[3] for ss in flux_info][13:]
        flux430  = [ss[4] for ss in flux_info][13:]
        flux820  = [ss[5] for ss in flux_info][13:]
        flux1380 = [ss[6] for ss in flux_info][13:]
        flux1500 = [ss[7] for ss in flux_info][13:]
        flux2000 = [ss[8] for ss in flux_info][13:]
        names = [ss[0] for ss in flux_info][13:]
        fig2 = plt.figure()
        gs2 = gridspec.GridSpec(8, 16)
    if i==2:
        flux57   = [[ss[1] for ss in flux_info][7]]
        flux327  = [[ss[2] for ss in flux_info][7]]
        flux350  = [[ss[3] for ss in flux_info][7]]
        flux430  = [[ss[4] for ss in flux_info][7]]
        flux820  = [[ss[5] for ss in flux_info][7]]
        flux1380 = [[ss[6] for ss in flux_info][7]]
        names = [[ss[0] for ss in flux_info][7]]
        fig3 = plt.figure()
        
    periods = []
    dms = []

    for n in names:
        par_path = f"{DATA_PATH}{n}_fiore+22.par"
        with open(par_path, 'r') as infile:
            for l in infile.readlines():
                if l.startswith("F0"):
                    f0 = float(l.split()[1])*u.Hz
                    p0_s = 1./f0
                    if p0_s < 10.*u.ms:
                        p0_ms = round(p0_s.to(u.ms).value,1)
                        p0_label = [p0_ms,"ms"]
                    elif p0_s > 999.99*u.ms:
                        p0_label = [round(p0_s.value,2),"s"]
                    else:
                        p0_s = 1./f0
                        p0_ms = int(p0_s.to(u.ms).value)
                        p0_label = [p0_ms,"ms"]
                    periods.append(p0_label)
                elif l.startswith("DM "):
                    dm = float(l.split()[1])
                    if dm<10.:
                        dm = round(dm,1)
                    else:
                        dm = int(dm)
                    dms.append(dm)

    per_u = [f"{pp[0]} {pp[1]}" for pp in periods]
    dm_u = [f"{dd} pc/cm$^3$" for dd in dms]

    nbin = 128
    midbin = nbin/2
    both = 0
    count = 0
    phase = np.arange(nbin)/float(nbin)

    print(names)
    for nn,pp,dd,s57,s327,s350,s430,s820,s1380,s1500,s2000 in zip(names,per_u,dm_u,flux57,flux327,flux350,flux430,flux820,flux1380,flux1500,flux2000):

        profile_fname_57  = f"{DATA_PATH}{nn}_57.15MHz_fiore+22.profile"
        profile_fname_327 = f"{DATA_PATH}{nn}_327MHz_fiore+22.profile"
        profile_fname_350 = f"{DATA_PATH}{nn}_350MHz_fiore+22.profile"
        profile_fname_430 = f"{DATA_PATH}{nn}_430MHz_fiore+22.profile"
        profile_fname_820 = f"{DATA_PATH}{nn}_820MHz_fiore+22.profile"
        profile_fname_1380 = f"{DATA_PATH}{nn}_1380MHz_fiore+22.profile"
        profile_fname_1500 = f"{DATA_PATH}{nn}_1500MHz_fiore+22.profile"
        profile_fname_2000 = f"{DATA_PATH}{nn}_2000MHz_fiore+22.profile"

        print(nn)
        row = int(np.floor(count/4))
        col = count % 4

        r_lo, r_hi, c_lo, c_hi = (4*row,4*(row+1), 4*col,4*(col+1))
        print(row, col)
        if i==0:
            ax1 = plt.subplot(gs1[r_lo:r_hi, c_lo:c_hi])
        elif i==1:
            ax2 = plt.subplot(gs2[r_lo:r_hi, c_lo:c_hi])
        else:
            ax3 = plt.subplot()

        shift  = 0.0
        offset = 0.0
        rotate = 0
        if nn=="J0636+5128":
            dshift = 0.2
        elif nn=="J1327+3423":
            shift = 0.0
            dshift = 0.15
            offset = -0.05
        elif nn=="J1434+7257":
            dshift = 0.25
        elif nn=="J1816+4510":
            dshift = 0.25
            rotate = 30-64
        elif nn=="J2210+5712":
            offset = 0.3
        else:
            dshift = 0.4

        if profile_fname_57 in all_files:
            y = prof_function(profile_fname_57)
            plt.plot(phase,np.roll(y,rotate)+offset,c='black',label="57 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_327 in all_files:
            y = prof_function(profile_fname_327)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='purple',label="327 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_350 in all_files:
            y = prof_function(profile_fname_350)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='r',label="350 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_430 in all_files:
            y = prof_function(profile_fname_430)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='pink',label="430 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_820 in all_files:
            y = prof_function(profile_fname_820)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='b',label="820 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_1380 in all_files:
            y = prof_function(profile_fname_1380)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='g',label="1380 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_1500 in all_files:
            y = prof_function(profile_fname_1500)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='g',label="1500 MHz" if count == 0 else "")
            shift += dshift
        if profile_fname_2000 in all_files:
            y = prof_function(profile_fname_2000)
            plt.plot(phase,np.roll(y,rotate)+shift+offset,c='cyan',label="2000 MHz" if count == 0 else "")
            shift += dshift

    #     # Center max profile value
    #     max_bin = np.argmax(x)+15
    #     int_shift_by = int(midbin-max_bin)
    #     x_r = np.roll(x,int_shift_by)
    #     if both: y_r = np.roll(y,int_shift_by)

    #     if nn == "J2022+2534":
    #         x_r = np.roll(x,int_shift_by-10)
    #         y_r = np.roll(y,int_shift_by+17)

    #     plt.plot(phase,x_r,c='r',label="350 MHz" if count == 0 else "")
    #     if both: plt.plot(phase,y_r+0.5,c='b',label="820 MHz" if count == 0 else "")
        plt.ylim([-0.1,1.25])
        count += 1

        
    #     plt.text(0.0,0.3,s350,horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='r')
    #     plt.text(0.0,0.15,"mJy",horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='r')
    #     if both: plt.text(0.0,0.88,s820,horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='b')
    #     if both: plt.text(0.0,0.73,"mJy",horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='b')

        if i==0:
            plt.text(0.08,1.08,f"PSR {nn.replace('-','$-$')}",horizontalalignment='left',verticalalignment='bottom',fontsize=8)
            plt.text(1.,0.95,dd,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            plt.text(0.95,0.8,pp,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            ax1.set_xticks([])
            ax1.set_yticks([])
        elif i==1:
            plt.text(0.05,1.1,f"PSR {nn.replace('-','$-$')}",horizontalalignment='left',verticalalignment='bottom',fontsize=9)
            plt.text(1.,0.95,dd,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            plt.text(0.95,0.8,pp,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            ax2.set_xticks([])
            ax2.set_yticks([])
        else:
            plt.text(0.05,1.1,f"PSR {nn.replace('-','$-$')}",horizontalalignment='left',verticalalignment='bottom',fontsize=9)
            plt.text(1.,0.95,dd,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            plt.text(0.95,0.8,pp,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
            ax3.set_xticks([])
            ax3.set_yticks([])
        
#plt.show()
fig1.savefig('profiles_1.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
fig2.savefig('profiles_2.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
fig3.savefig('J1327_profile.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
