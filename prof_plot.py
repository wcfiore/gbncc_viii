import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import glob

def prof_function(file):

    x,y,z,profile = np.loadtxt(file,dtype='float',unpack=True)

    # Roughly identify off-pulse region
    init_median = np.median(profile)
    off_bins = np.where(profile < init_median)
    off_mean = np.mean(profile[off_bins])

    # Scale profile so max = 0.5 and off-pulse noise is centered on 0.
    profile -= off_mean
    profile /= 2 * np.max(profile) 

    return profile

fig = plt.figure()
gs = gridspec.GridSpec(12, 16)


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
DATA_PATH = "data/clutter/old_profiles/" # using (nearly identical) old profiles since new ones are not aligned
all_files = glob.glob(f"{DATA_PATH}*.profile")
flux_info = np.loadtxt(f"{DATA_PATH}flux.info",dtype="str")
flux350 = [ss[1] for ss in flux_info]
flux820 = [ss[2] for ss in flux_info]
names = [ss[0] for ss in flux_info]

periods = ['64','3.1','83','24','7.8','1.9','2.9','132','225','41','2.6','3.3']
dms = ['54','21','17','4.8','40','16','29','82','26','31','54','24']
per_u = [f"{pp} ms" for pp in periods]
dm_u = [f"{dd} pc/cm$^3$" for dd in dms]

nbin = 128
midbin = nbin/2
both = 0
count = 0
phase = np.arange(nbin)/float(nbin)

print(names)
for nn,pp,dd,s350,s820 in zip(names,per_u,dm_u,flux350,flux820):

    profile_fname_350 = f"{DATA_PATH}{nn}_350MHz_swiggum+22.profile"
    profile_fname_820 = f"{DATA_PATH}{nn}_820MHz_swiggum+22.profile"

    print(nn)
    row = int(np.floor(count/4))
    col = count % 4

    r_lo, r_hi, c_lo, c_hi = (4*row,4*(row+1), 4*col,4*(col+1))
    print(row, col)
    ax = plt.subplot(gs[r_lo:r_hi, c_lo:c_hi])

    x = prof_function(profile_fname_350)

    if profile_fname_820 in all_files:
        both = 1
    else:
        both = 0

    if both: y = prof_function(profile_fname_820)

    # Center max profile value
    max_bin = np.argmax(x)+15
    int_shift_by = int(midbin-max_bin)
    x_r = np.roll(x,int_shift_by)
    if both: y_r = np.roll(y,int_shift_by)

    if nn == "J2022+2534":
        x_r = np.roll(x,int_shift_by-10)
        y_r = np.roll(y,int_shift_by+17)

    plt.plot(phase,x_r,c='r',label="350 MHz" if count == 0 else "")
    if both: plt.plot(phase,y_r+0.5,c='b',label="820 MHz" if count == 0 else "")
    plt.ylim([-0.1,1.25])
    count += 1

    plt.text(0.05,1.1,f"PSR {nn.replace('-','$-$')}",horizontalalignment='left',verticalalignment='bottom',fontsize=9)
    plt.text(0.97,0.9,dd,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
    plt.text(0.92,0.75,pp,horizontalalignment='right',verticalalignment='bottom',fontsize=6)
    plt.text(0.0,0.3,s350,horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='r')
    plt.text(0.0,0.15,"mJy",horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='r')
    if both: plt.text(0.0,0.88,s820,horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='b')
    if both: plt.text(0.0,0.73,"mJy",horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='b')

    ax.set_xticks([])
    ax.set_yticks([])

fig.savefig('profiles.pdf',format='pdf',bbox_inches='tight',pad_inches=0.25)
