"""
Plot 2.3m ascii spectra from a specified folder
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

try:
    root = sys.argv[1]
except:
    root = '/Users/marusa/observing/23m/data/ascii/'

filenamesB = []
for f in os.listdir(root):
    if f.endswith("_b.dat"):
        filename = os.path.join(root, f)
        filenamesB.append(filename)

filenamesR = []
for f in os.listdir(root):
    if f.endswith("_r.dat"):
        filename = os.path.join(root, f)
        filenamesR.append(filename)

filenamesB = sorted(filenamesB)
filenamesR = sorted(filenamesR)

print len(filenamesB)
print len(filenamesR)

##############################
###### RED ###################
##############################
li_min=6705
li_max=6711

cont_min = 6690
cont_max = 6720

fig=plt.figure()
ax=fig.add_subplot(111)
offset=1
postext=6712
#~ clr = cycle(['k', 'b'])
#~ c=clr.next()

# Read and normalise spectra
for i, filename in enumerate(filenamesR):
    if '20190228' in f:
        continue

    d=np.loadtxt(filename)
    d0=d

    #~ # Quick quasi normalization
    mask = (d[:,0]>cont_min) & (d[:,0]<cont_max)
    d=d[mask]
    
    # exclude lithium
    mask_cont = np.logical_or(d[:,0]<li_min, d[:,0]>li_max)

    median = np.nanmedian(d[mask_cont,1])
    #~ d[:,1] = d[:,1]/median

    w = d0[:,0]
    f = d0[:,1]/median

    
    ax.plot(w, f+offset*i, c='k', linewidth=0.4, alpha=1)
    
    #~ ax.annotate('%.2f'%ew, xy=(postext, 1+offset*j), xytext=(postext, 1+offset*j))
    
    i+=1



#~ oprev=filename0.split('-')[1].split('.')[0]


#~ ax2.set_xlim(6520, 6600)
#~ ax.set_ylim(-10, 160)
#~ ax.set_ylim(0.7, 2.5)
#~ ax.set_ylim(0, 13)
#~ ax.set_xlim(6520, 6600)
ax.axvline(x=6708, color='r', linewidth=0.3)
ax.axvline(x=li_min, color='k', linewidth=0.3)
ax.axvline(x=li_max, color='k', linewidth=0.3)


ax.axvline(x=6563, color='k', linewidth=0.3)
ax.axvline(x=6563-2.7, color='k', linewidth=0.3)
ax.axvline(x=6563+2.7, color='k', linewidth=0.3)

plt.tight_layout()

li_min=6705
li_max=6711

cont_min = 6690
cont_max = 6720



##############################
###### BLUE ##################
##############################
fig=plt.figure()
ax=fig.add_subplot(111)
offset=1
postext=6712


# Read and normalise spectra
for i, filename in enumerate(filenamesB):
    if '20190228' in f:
        continue

    d=np.loadtxt(filename)
    d0=d

    #~ # Quick quasi normalization
    #~ mask = (d[:,0]>cont_min) & (d[:,0]<cont_max)
    #~ d=d[mask]
    
    # exclude lithium
    #~ mask_cont = np.logical_or(d[:,0]<li_min, d[:,0]>li_max)

    #~ median = np.nanmedian(d[mask_cont,1])
    #~ d[:,1] = d[:,1]/median

    w = d0[:,0]
    f = d0[:,1]#/median

    
    ax.plot(w, f+offset*i, c='k', linewidth=0.4, alpha=1)
    
    #~ ax.annotate('%.2f'%ew, xy=(postext, 1+offset*j), xytext=(postext, 1+offset*j))
    
    i+=1

plt.tight_layout()


plt.show()
