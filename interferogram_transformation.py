## Imported modules ##
import numpy as np
import matplotlib.pyplot as plt
## Path where is located interferograms ##
spt_path='/home/STG05_TCCON/ggg2020/spt/'

## Name of files that contain data ##
filename=spt_path+'zpa20040721saaaaa.043'
filename=spt_path+'zoc20110208lfddaa.001'

## Manipulation of dataset## 
spcmatrix = np.genfromtxt(filename,skip_header=3)
print (spcmatrix)
fig = plt.figure()
Tm = spcmatrix[:,1]
Tc = spcmatrix[:,2]
r = Tm - Tc
## Interferogram plotting ##
ax = fig.add_subplot(1,1,1)
ax.set_xlabel("Wavenumber (cm-1)")
ax.set_ylabel("Transmittance")
plt.plot(spcmatrix[:,0],spcmatrix[:,1],label="measurement" )
plt.plot(spcmatrix[:,0],spcmatrix[:,2],label="calculation")
#plt.plot(spcmatrix[:,0], r ,label="residual")
ax.legend(loc="best")
plt.show()
