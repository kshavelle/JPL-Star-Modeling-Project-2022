#Package libraries
from radmc3dPy.image import *
from matplotlib import cm

#Plotting the image
makeImage(npix=200.0,incl=60.0,phi=30.0,wav=1.3e3,sizeau=200.0)
#npix = Specifies the number of pixels used to make the image, this is adjustable. Here, I use 200, meaning an image of 200x200 pixels will be made. You can also specify the x- and y- direction number of pixels separately, using "npixx" and "npixy."
#incl = This specifies the inclination angel the object is viewed at in degrees, this is adjustable.
#phi = Gives the remaining angel in degrees, related to inclination. For example, if incl=90 and phi=0, this means that the observer is located at infinity toward the negative y axis.
#wav = Gives the wavelength the object is imaged at.
#sizeau = Provides the image size in model space in units of AU (=1.496E13 cm).
im_mm = readImage()
plt.figure()
plotImage(im_mm,au=True,log=True,maxlog=3,bunit='jy/pixel',dpc=dpc,cmap='magma')
#dpc = Distance of observer in units of parsec.
#bunit = Units of the image.
#cmap = Matplotlib colormap.


