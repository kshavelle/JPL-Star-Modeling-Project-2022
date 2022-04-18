#Package libraries
from radmc3dPy.analyze import *
import matplotlib.pyplot as plt

#Creating the SED
tag = 'Disk + Bipoutflow with Envelope '
s = readSpectrum(fname = 'spectrum.out') #column 0: wavelength in microns; column 1: Flux in cgs.
distance = 710.0 #in pc. The spectrum.out file is still normalized to a distance of 1 pc (see radmc3d docs)
F_nu = s[:,1.0] * distance**-2.0 * 1.0e23 #to Jy at the set distance
nu = 3e8 * s[:,0]**-1.0 * 1.0e6 * 1.0e-9 #microns to GHz
plt.plot(nu, F_nu)
plt.title('%s - distance: %d pc'%(tag,distance))
plt.xlabel('Frequency [GHz]'); plt.ylabel('Flux density [Jy]')
plt.xscale('log'); plt.yscale('log')
plt.savefig('sed_'+tag+'.png')
plt.show()
