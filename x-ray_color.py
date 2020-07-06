#""" uses  
#         
#	 
#"""
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

# read image into numpy-array ima
imarray = mpimg.imread(outputname+'.tif')
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

