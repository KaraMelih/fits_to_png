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


args = parser.parse_args()

# give meaningful variablenames to the command line
# arguments:
fits = args.fits[0]
outputname = args.output_file[0]
conf = args.conf[0]
cmap = args.cmap[0]
ip = args.interpolation[0]
sm = args.smooth[0]
cr = args.crop[0]
br = float(args.brightness[0])
cmap_check = args.cmap_check[0]
smooth_check = args.smooth_check[0]

if '.' in outputname:
    outputname = outputname.split('.')[0]

# set the coordinate system
hdu = FITS.open(fits)
wcs = WCS(hdu[0].header)

try:
    sp.call(["stiff", f"{fits}","-OUTFILE_NAME", f"{outputname}.tif","-c",f"{conf}"]);
except:
    print('Configuration file not found. Using default.')
    sp.call(["stiff", f"{fits}","-OUTFILE_NAME", f"{outputname}.tif"]);

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
        cmap = plt.get_cmap(cm)
        colored_img = cmap(imarray)
        ax.imshow(colored_img[::binn,::binn]*br, interpolation=ip)
        ax.set_title(cm)
        ax.axis('off')
    figure.tight_layout() 
    plt.savefig(outputname+'_alternative_cmaps.png')
    
def smooth_alternatives(imarray,br=br,ip=ip, binn=1, cm=cmap):
    '''
    creates a 3x3 plots for different smoothing alternatives
    '''
    from scipy.ndimage import gaussian_filter
    from scipy.ndimage import median_filter
    
    imarray = imarray[::binn,::binn]
    
    cmap = plt.get_cmap(cm)
    colored_img = cmap(imarray)
    
    median_filtered = []
    for window_size in [1,3,5]:
        ims = []
        for d in range(3):
            im_conv_d = median_filter(colored_img[:,:,d], size=(d**2+window_size,d**2+window_size), mode='reflect')
            ims.append(im_conv_d)

        im_conv = np.stack(ims, axis=2) #.astype("uint8")
        median_filtered.append(im_conv)

    smoothed = []
    for sm in [1.5,3,5]:
        im = gaussian_filter(imarray, sigma=sm)
        im = cmap(im)
        smoothed.append(im)
        
    alternatives = median_filtered + smoothed
    titles = ['median_filter:1','median_filter:3','median_filter:5','gaussian_filter:1.5','gaussian_filter:3','gaussian_filter:5',]
    figure, [[a11,a12,a13],[a21,a22,a23]] = plt.subplots(nrows=2, ncols=3,figsize=(15,10))
    for i,(sm,ax) in enumerate(zip(alternatives,[a11,a12,a13,a21,a22,a23])):
        ax.imshow(sm, interpolation=ip)
        ax.set_title(titles[i])
        ax.axis('off')
    figure.tight_layout() 
    plt.savefig(outputname+'_alternative_smoothing.png')

# read image into numpy-array ima
imarray = mpimg.imread(outputname+'.tif')

if (cmap_check != 'N') or (smooth_check != 'N'):
    if cmap_check != 'N':
        if sm != 'None':
            from scipy.ndimage import gaussian_filter
            imarray = gaussian_filter(imarray, sigma=float(sm))
        cmap_alternatives(imarray, binn=10) # binning so its faster
    if smooth_check != 'N':
        smooth_alternatives(imarray,binn=10)
    sys.exit()

if sm != 'None':
    from scipy.ndimage import gaussian_filter
    imarray = gaussian_filter(imarray, sigma=float(sm))

# color the image (i.e. add rgba axis)
cmap = plt.get_cmap(cmap)
# print(imarray.shape)
colored_img = cmap(imarray)
# print(colored_img.shape)
    
#  draw figure
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, projection=wcs)
ax.imshow(np.flipud(colored_img*br), interpolation=ip)

if cr == 'True':
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
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

