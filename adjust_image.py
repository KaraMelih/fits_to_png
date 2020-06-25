#!/usr/bin/env python3

# HISTORY
# 08.02.2020 Created
# 20.03.2020
# I added the new NED redshifts (only for galaxies)
# 
# 20.06.2020
# I made some additions, added binning.
# Now it is compatible with create_subimages.py 
# Still, it needs to be adjusted such that it can overlay ellipse regions.

# standard libraries
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as au
import astropy.coordinates as ac
import astropy.wcs.utils as awu
from PIL import Image
import argparse
import sys
import os
import matplotlib.patheffects as path_effects
import warnings
warnings.filterwarnings("ignore")

# Handle command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
DESCRIPTION:
The program takes an image file, a fits file for coordinate information
and ds9 contours file as input
creates high quality contours overlaid images.

Required inputs (-img) image file, a fits file (-f), (-con) contours to overlay 
(-o) output figure name.

EXAMPLES:
- adjust_image.py -img stiff.tif -f r_cut.fits -con mycontour.con -o finalimg.png
  Takes stiff.tif image plots according to wcs of r_cut.fits file,
  overlays mycontour.con contours and creates finalimg.png

- adjust_image.py -img stiff.tif -f r_cut.fits -reg ds9regions.reg -o myimg.tif
  Takes stiff.tif image plots according to wcs of r_cut.fits file,
  overlays ds9regions.reg regions and creates myimg.tif

AUTHOR:
Melih Kara (karamel@astro.uni-bonn)


"""
)
parser.add_argument('-img', '--input_img', nargs=1,
                    required=True, help='input image name')
parser.add_argument('-f', '--fits_img', nargs=1,
                    required=True, help='input fits image')
parser.add_argument('-con', '--contour_file', nargs=1, required=False, default=['None'],
                    type=str, help='contour file')
parser.add_argument('-o', '--output_file', nargs=1,
                    required=True, help='output file ')
parser.add_argument('-col', '--color', nargs=1, default=['white'],
                    type=str, help='contour color (OPTIONAL) \n ex: -col blue')
parser.add_argument('-lw', '--line_width', nargs=1, default=[1],
                    type=float, help='contour line width (OPTIONAL) \n ex: -lw 2')
parser.add_argument('-reg', '--region_file', nargs=1, default=['None'],
                    type=str, help='region file (OPTIONAL) \n ex: -reg ds9.reg')
parser.add_argument('-reg2', '--region_file2', nargs=1, default=['None'],
                    type=str, help='region file (OPTIONAL) \n ex: -reg2 ds9_new.reg')
parser.add_argument('-reg3', '--region_file3', nargs=1, default=['None'],
                    type=str, help='region file (OPTIONAL) \n ex: -reg3 ds9_new3.reg')
parser.add_argument('-dpi', '--resolution', nargs=1, default=[300],
                    type=int, help='resolution (OPTIONAL) \n ex: -dpi 1200')
parser.add_argument('-magic', '--magic', nargs=1, default=[0],
                    type=int, help='Overlay known-z objects, only works for A3391/95 system (OPTIONAL) \n ex: -magic 1')
parser.add_argument('-bin', '--binning', nargs=1, default=[1],
                    type=int, help='Choose a binning, default is 1. For large data resolution \
                                    can be decreased by binning the data (OPTIONAL) \n ex: -bin 5')

args = parser.parse_args()

# give meaningful variablenames to the command line
# arguments:
input_image = args.input_img[0]
input_fits = args.fits_img[0]
output_fig = args.output_file[0]
input_cont = args.contour_file[0]
color = args.color[0]
lw = args.line_width[0]
regionfile = args.region_file[0]
regionfile2 = args.region_file2[0]
regionfile3 = args.region_file3[0]
dpi = args.resolution[0]
magic = args.magic[0]
binning = args.binning[0]

# get WCS from FITS image:
hdu = fits.open(input_fits)
wcs = WCS(hdu[0].header)

# Create a figure
# projection is with respect to the fits file coordinates
fig = plt.figure(figsize=(10,10))
plt.style.use('seaborn-colorblind')
ax = plt.subplot(projection=wcs)
ax.set_xlim(-0.5, hdu[0].data.shape[1]/binning + 0.5)
ax.set_ylim(-0.5, hdu[0].data.shape[0]/binning + 0.5)
ra = ax.coords[0]
dec = ax.coords[1]

if binning == 1:
    ra.set_major_formatter('d.dd')
    dec.set_major_formatter('d.dd')

ax.set_xlabel('RA',fontsize=15, rotation=270)
ax.set_ylabel('DEC',fontsize=15)
ax.tick_params(axis='both', which='both', labelsize=20)

if magic ==1 :
# Display known redshifts
    posx,posy,reds, typs = np.loadtxt('/vol/aibn1058/data1/karamel/Thesis/DECam_all/extended_srcs_MRC/NED_data_Kara.txt', unpack=True)
    posx,posy,reds = np.loadtxt('/vol/aibn1058/data1/karamel/Thesis/DECam_all/extended_srcs_MRC/NED_data_new20-03-2020.txt', unpack=True)
    xp = posx * au.deg
    yp = posy * au.deg
    cp = ac.SkyCoord(ra=xp, dec=yp)
    ra_pixned, dec_pixned = awu.skycoord_to_pixel(cp, wcs)
##print(posx[0],posy[0],reds[0])
    for i in range(len(ra_pixned)):
        text = ax.text(ra_pixned[i],dec_pixned[i], str(reds[i]), color='white', ha='center', va='center', size=15, clip_on=True)
        text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])


# Display the image
img = np.array(Image.open(input_image))
# img.setflags(write=1)
img[:,:,0] = np.flipud(img[:,:,0])
img[:,:,1] = np.flipud(img[:,:,1])
img[:,:,2] = np.flipud(img[:,:,2])

ax.imshow(img, origin='lower')

# Overlay the contours
# Open contour files and select the levels
if input_cont != 'None':
#    import analysis_funcs as af
#    af.set_contour(input_cont,ax, wcs=wcs, alpha=1, lw=25, color="white")
    cont_orig = open(input_cont)
    lines = cont_orig.readlines()
    lines=np.array(lines)
    cont_levels = np.array(np.where(lines=='\n'))[0]
    # make a temporary file
    lines_new = lines.copy()
    lines_new[cont_levels] = ' 0 0 \n'
    temp_con = open('./tmp_con.con','w+')
    temp_con.writelines(lines_new)
    temp_con.close()
    cont = np.loadtxt('./tmp_con.con')
    for i in range(2,len(cont_levels)-1):    
        x = cont[cont_levels[i]+1:cont_levels[i+1]][:,0] * au.deg
        y = cont[cont_levels[i]+1:cont_levels[i+1]][:,1] * au.deg
        c = ac.SkyCoord(ra=x, dec=y)

        ra_pix, dec_pix = awu.skycoord_to_pixel(c, wcs)
        ax.plot(ra_pix/binning, dec_pix/binning, '-', color=color, lw=lw)
    os.remove('./tmp_con.con')

#pixel = 7.311111523045E-05 # degrees
if regionfile != 'None':
    regfile = open(regionfile)
    reglines = regfile.readlines()
    reglines = np.array(reglines)
    reglines = reglines[3:]
    RAreg = np.array([float(i.split('(')[1].split(',')[0]) for i in reglines])
    DECreg = -1*np.array([float(i.split('-')[1].split(',')[0]) for i in reglines])
    rad = np.array([float(i.split(',')[2].split('"')[0]) for i in reglines])

    x1 = RAreg * au.deg
    y1 = DECreg * au.deg
    c1 = ac.SkyCoord(ra=x1, dec=y1)
    ra_pix1, dec_pix1 = awu.skycoord_to_pixel(c1, wcs)
    y2 = (DECreg+rad/3600) * au.deg
    c2 = ac.SkyCoord(ra=x1, dec=y2)
    _, dummy = awu.skycoord_to_pixel(c2,wcs)
    rad_pix = abs(dummy-dec_pix1)
    for i in range(len(ra_pix1)):
        circle = plt.Circle((ra_pix1[i]/binning, dec_pix1[i]/binning), rad_pix[i]/binning, color='white', fill=False, lw=1)
        ax.add_artist(circle)


if regionfile2 != 'None':
    regfile = open(regionfile2)
    reglines = regfile.readlines()
    reglines = np.array(reglines)
    reglines = reglines[3:]
    RAreg = np.array([float(i.split('(')[1].split(',')[0]) for i in reglines])
    DECreg = -1*np.array([float(i.split('-')[1].split(',')[0]) for i in reglines])
    rad = np.array([float(i.split(',')[2].split('"')[0]) for i in reglines])
    x1 = RAreg * au.deg
    y1 = DECreg * au.deg
    c1 = ac.SkyCoord(ra=x1, dec=y1)
    ra_pix1, dec_pix1 = awu.skycoord_to_pixel(c1, wcs)

    y2 = (DECreg+rad/3600) * au.deg
    c2 = ac.SkyCoord(ra=x1, dec=y2)
    _, dummy = awu.skycoord_to_pixel(c2,wcs)
    for i in range(len(ra_pix1)):
        circle = plt.Circle((ra_pix1[i], dec_pix1[i]), rad_pix[i], color='purple', fill=False, lw=1)#, transform=ax.get_transform('fk5')
        ax.add_artist(circle)

if regionfile3 != 'None':
    regfile = open(regionfile3)
    reglines = regfile.readlines()
    reglines = np.array(reglines)
    reglines = reglines[3:]
    RAreg = np.array([float(i.split('(')[1].split(',')[0]) for i in reglines])
    DECreg = -1*np.array([float(i.split('-')[1].split(',')[0]) for i in reglines])
    rad = np.array([float(i.split(',')[2].split('"')[0]) for i in reglines])
    x1 = RAreg * au.deg
    y1 = DECreg * au.deg
    c1 = ac.SkyCoord(ra=x1, dec=y1)
    ra_pix1, dec_pix1 = awu.skycoord_to_pixel(c1, wcs)

    y2 = (DECreg+rad/3600) * au.deg
    c2 = ac.SkyCoord(ra=x1, dec=y2)
    _, dummy = awu.skycoord_to_pixel(c2,wcs)
    for i in range(len(ra_pix1)):
        circle = plt.Circle((ra_pix1[i], dec_pix1[i]), rad_pix[i], color='orange', fill=False, lw=1, ls='--') #, transform=ax.get_transform('fk5')
        ax.add_artist(circle)

#plt.savefig(output_fig, bbox_inches='tight', pad_inches=0.2, format='tiff', dpi=dpi)
plt.savefig(output_fig, bbox_inches='tight', format='tiff', pad_inches=0.2, dpi=dpi)
plt.close()
plt.cla()
plt.clf()



#ra.set_major_formatter('x.xxx')  # decimal, non-angle coordinates,
                                  # 3 decimal places

#ra.set_major_formatter('d')  # decimal degrees, no decimal places

#ra.set_major_formatter('d.ddddd')  # decimal degrees, 5 decimal places

#ra.set_major_formatter('dd:mm')  # sexagesimal, 1' precision
#ra.set_major_formatter('dd:mm:ss')  # sexagesimal, 1" precision
#ra.set_major_formatter('dd:mm:ss.ss')  # sexagesimal, 0.01" precision

#ra.set_major_formatter('hh:mm')  # sexagesimal (hours)
#ra.set_major_formatter('hh:mm:ss')  # sexagesimal (hours)
#ra.set_major_formatter('hh:mm:ss.ss')  # sexagesimal (hours)










