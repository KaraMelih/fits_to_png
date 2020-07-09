# HISTORY  
#         
#	 
# Currently, median_filter in the general case is not implemented
# 

# Libraries
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.io import fits as FITS
from astropy.wcs import WCS
import matplotlib.image as mpimg
import subprocess as sp
import argparse
import sys
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter

import warnings
warnings.filterwarnings("ignore")

# set plotting params
matplotlib.rcParams['font.size']=18.0
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.labelsize'] = 16.0
matplotlib.rcParams['ytick.labelsize'] = 16.0

# Handle command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
DESCRIPTION:
    Script to create nicely colored X-ray images.

Required inputs (-f) image fits file, (-o) output name
Optional [[-sm smoothing radius], [-cmap colormap], [-ip interpolation method], [-cr crop image axis]]

EXAMPLES: 

 - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -camp jet
 o Creates a grey colored stiff image named mycolored_image.tif, colors it with jet color map and saves as mycolored_image.png 
 
 - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -camp jet -br 1.7
 o Creates the same image, increases brightness
 
 - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -camp jet -br 1.7 -sm 3
 o Creates the same image, increases brightness, gaussian smooths with a 3 sigma kernel
 
  - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -camp jet -br 1.7 -sm 3 -cr True
 o Creates the same image, increases brightness, gaussian smooths with a 3 sigma kernel, removes the axis names, labels, ticks. 
 
  - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -camp jet -sc Y 
 o creates a mycolored_image_smoothing_alternatives.png with a several different smoothing alternatives for the given colormap. Exits after that

  - python3 x-ray_color.py -f myimgae.fits -o mycolored_image -sm 3 -cc Y
 o creates a mycolored_image_cmap_alternatives.png with a several different colormap alternatives for the given smothing. Exists after that. 

AUTHOR:
Melih Kara (karamel@astro.uni-bonn)


"""
)
parser.add_argument('-f', '--fits', nargs=1,
                    required=True, help='input image name')
parser.add_argument('-o', '--output_file', nargs=1,
                    required=True, help='output file ')
parser.add_argument('-cmap', '--cmap', nargs=1,
                    required=False,default=['cubehelix'], help='color map, deafult: cubehelix')
parser.add_argument('-ip', '--interpolation', nargs=1,
                    required=False,default=['bicubic'], help='interpolation, deafult: bicubic')
parser.add_argument('-sm', '--smooth', nargs=1,
                    required=False,default=['None'], help='Gaussian smoothing eg: sm 3 , default is None')
parser.add_argument('-sm_method', '--smooth_method', nargs=1,
                    required=False, default=['G'], help='Smoothing method Gaussian(G) or Median(M) default is Gaussian \n ex: -sm_method G ')
parser.add_argument('-cr', '--crop', nargs=1,
                    required=False,default=['False'], help='Crop the axis, labels, etc eg: cr True , default is False')
parser.add_argument('-br', '--brightness', nargs=1,
                    required=False, default=[1], help='Increase/Decrease brightness \nEX: br 1.6 , default is 1')
parser.add_argument('-conf', '--conf', nargs=1,
                    required=False,default=['x_ray.conf'], help='')
parser.add_argument('-cc', '--cmap_check', nargs=1,
                    required=False,default=['N'], help='bring several cmap alternatives. syntax: -cc Y')
parser.add_argument('-sc', '--smooth_check', nargs=1,
                    required=False,default=['N'], help='bring several smoothing alternatives. syntax: -sc Y')
parser.add_argument('-vmin', '--vmin', nargs=1,
                    required=False, default=[0.075], help='min value in logarithmic plot ex: -vmin 0.1')


args = parser.parse_args()

# give meaningful variablenames to the command line
# arguments:
fits = args.fits[0]
outputname = args.output_file[0]
conf = args.conf[0]
cmap = args.cmap[0]
ip = args.interpolation[0]
sm = args.smooth[0]
sm_method = args.smooth_method[0].lower()
cr = args.crop[0]
br = float(args.brightness[0])
cmap_check = args.cmap_check[0]
smooth_check = args.smooth_check[0]
vmin = args.vmin[0]

if sm_method != 'g' or sm_method != 'm':
    print('Undefined smooth method! Only "G" or "M" allowed\n Using Gaussian')
    sm_method = 'g'

if '.' in outputname:
    outputname = outputname.split('.')[0]

# set the coordinate system
hdu = FITS.open(fits)
wcs = WCS(hdu[0].header)         # world-coordinate-system
img = hdu[0].data                # image data

# define some initial-check functions
def cmap_alternatives(imarray,br=br,ip=ip, binn=1):
    '''
    creates a 3x3 plots for different color map alternatives
    '''
    cmap_alternatives = ['twilight_r','plasma','cubehelix',\
                        'jet','gist_earth','gist_stern',\
                        'gist_ncar','rainbow','CMRmap']    
    
    figure, [[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]] = plt.subplots(nrows=3, ncols=3,figsize=(15,15))
    for cm,ax in zip(cmap_alternatives,[a11,a12,a13,a21,a22,a23,a31,a32,a33]):
        ax.imshow(np.flipud(imarray[::binn,::binn])*br, interpolation=ip, cmap=cm, norm=LogNorm(vmin=vmin, vmax=np.max(imarray)))
        ax.set_title(cm)
        ax.axis('off')
    figure.tight_layout() 
    plt.savefig(outputname+'+alternative_cmaps.png')
    
def smooth_alternatives(imarray,br=br,ip=ip, binn=1, cm='cubehelix'):
    '''
    creates a 3x3 plots for different smoothing alternatives
    '''    
    imarray = imarray[::binn,::binn]
    
    median_filtered = []
    for window_size in [1,3,5]:
        im = median_filter(imarray, size=window_size)
        median_filtered.append(im)

    gaussian_filtered = []
    for sm in [1.5,3,5]:
        im = gaussian_filter(imarray, sigma=sm)
        gaussian_filtered.append(im)
        
    alternatives = median_filtered + gaussian_filtered
    titles = ['median_filter:1','median_filter:3','median_filter:5','gaussian_filter:1.5','gaussian_filter:3','gaussian_filter:5',]
    figure, [[a11,a12,a13],[a21,a22,a23]] = plt.subplots(nrows=2, ncols=3,figsize=(15,10))
    for i,(sm,ax) in enumerate(zip(alternatives,[a11,a12,a13,a21,a22,a23])):
        ax.imshow(np.flipud(sm), cmap=cm, interpolation=ip, norm=LogNorm(vmin=vmin, vmax=np.max(imarray)))
        ax.set_title(titles[i])
        ax.axis('off')
    figure.tight_layout() 
    plt.savefig(outputname+'+alternative_cmaps.png')


if (cmap_check != 'N') or (smooth_check != 'N'):
    if cmap_check != 'N':
        if sm != 'None':
            if sm_method == 'g'
                imarray = gaussian_filter(img, sigma=float(sm))
            else:
                imarray = median_filter(img, size=float(sm))
        cmap_alternatives(img, binn=10) # binning so its faster
    if smooth_check != 'N':
        smooth_alternatives(img,binn=10)
    sys.exit()

if sm != 'None':
    if sm_method == 'g'
        img = gaussian_filter(img, sigma=float(sm))
    else:
        img = median_filter(img, size=float(sm))


#  draw figure
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=wcs)
# show the image
ax.imshow(img*br, interpolation=ip, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=np.max(imarray)))

if cr == 'True':
    ax.tick_params(axis='x', which='both', bottom=False,  top=False, labelbottom=False)
    ax.tick_params(axis='y', which='both', bottom=False,  top=False, labelbottom=False) 
    plt.tight_layout() 
else:
    raax = ax.coords[0]
    decax = ax.coords[1]
    raax.set_major_formatter('d.d')
    decax.set_major_formatter('d.d')
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('DEC [deg]')

ax.grid(False)
plt.grid(b=None)

plt.savefig(outputname+'.png')

