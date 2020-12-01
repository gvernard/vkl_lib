import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import interpolate



b1_min = 4.7035368
b1_max = 8.2338171
b2_min = 4.0344134
b2_max = 7.1372285

Npix_b1 = 85
Npix_b2 = 50


x = np.linspace(b1_min,b1_max,Npix_b1)
y = np.linspace(b2_min,b2_max,Npix_b2)
xx, yy = np.meshgrid(x,y,sparse=True)
z  = 2*np.cos(xx*yy/3.0) + np.cos(5.0*xx)
zx = -2*yy*np.sin(xx*yy/3.0)/3.0 - 5*np.sin(5*xx)
zy = -2*np.sin(xx*yy/3.0)*xx/3.0
zxy = -2*np.sin(xx*yy/3.0)/3.0-2.0*np.cos(xx*yy/3.0)*xx*yy/9.0
zxx = -2*np.cos(xx*yy/3.0)*xx*yy/9.0 - 25*np.cos(5*xx)
zyy = -2*np.cos(xx*yy/3.0)*xx*xx/9.0


hdu = fits.PrimaryHDU(z)
hdu.writeto('python_z.fits',overwrite=True)

hdu = fits.PrimaryHDU(zx)
hdu.writeto('python_zx.fits',overwrite=True)

hdu = fits.PrimaryHDU(zx)
hdu.writeto('python_zy.fits',overwrite=True)

hdu = fits.PrimaryHDU(zxy)
hdu.writeto('python_zxy.fits',overwrite=True)

hdu = fits.PrimaryHDU(zxx)
hdu.writeto('python_zxx.fits',overwrite=True)

hdu = fits.PrimaryHDU(zyy)
hdu.writeto('python_zyy.fits',overwrite=True)


'''
fig, ax = plt.subplots(figsize=(10,8))
ax.pcolor(x,y,z)
#ax.set_aspect(1.0)
ax.set_xlim(b1_min,b1_max)
ax.set_ylim(b2_min,b2_max)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Function z')
plt.savefig('example_function.png')
'''

'''
fig, ax = plt.subplots(figsize=(10,8))
ax.pcolor(x,y,zx)
#ax.set_aspect(1.0)
ax.set_xlim(b1_min,b1_max)
ax.set_ylim(b2_min,b2_max)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Derivative z_x')
plt.savefig('example_function_dx.png')
'''
