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


N = 5
np.random.seed(1)
data = np.round(np.random.rand(N,N),1)
#print(data)

b1_min = -4.0
b1_max = 0.2
b2_min = -2.0
b2_max = 2.1

stepx_big = 0.84
stepy_big = 0.82
x = np.linspace(b1_min+stepx_big/2.0,b1_max-stepx_big/2.0,N)
y = np.linspace(b2_min+stepy_big/2.0,b2_max-stepy_big/2.0,N)
#print(x)
#print(y)
interp = 'cubic'
f = interpolate.interp2d(x,y,data,kind=interp)

stepx = 0.044
stepy = 0.042
xnew = np.linspace(b1_min+stepx_big/2.0,b1_max-stepx_big/2.0,30*N)
ynew = np.linspace(b2_min+stepy_big/2.0,b2_max-stepy_big/2.0,30*N)
znew = f(xnew,ynew)
#print(xnew)
#print(ynew)

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([0.125, 0.175, 0.75, 0.75])
plt.imshow(znew,extent=(xnew[0],xnew[-1],ynew[0],ynew[-1]),cmap='viridis',zorder=1)
plt.title(interp,weight='bold')
plt.xlim(xnew[0],xnew[-1])
plt.ylim(ynew[0],ynew[-1])

cax = fig.add_axes([0.125, 0.075, 0.75, 0.03])
cb = plt.colorbar(cax=cax,orientation='horizontal',ticks=np.linspace(0,1,6))
#cb.solids.set_edgecolor('face')

plt.scatter(xnew,ynew,marker='.',color='k',zorder=2,s=100)
plt.savefig("python_bilinear.png")


hdu = fits.PrimaryHDU(znew)
hdu.writeto('python_interp.fits',overwrite=True)
