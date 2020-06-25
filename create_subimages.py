#""" uses dfits_theli 
#         makesubimage
#	 stiff 
#"""
# assuming the reg files are circles with radii is written as inc arcsec with " at the end.

# standard libraries
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as au
import astropy.coordinates as ac
import astropy.wcs.utils as awu
import argparse
import sys
import os
import subprocess as sp
import warnings
warnings.filterwarnings("ignore")

# Handle command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
DESCRIPTION:
    Script to create and arrange sub-images from the given fits files
    Takes, 3-band fits files (optional) region file 
    If no region file is given, produces an image which is same size of the input
    For a single point region file, it will create an image around that region and place the region
    For multiple points region file it will create an image of the sizes of the input and overlay all the regions
    Additional region files can be feed, however all regions must have cirlce(deg,deg,arcsec") convention.
    CAVEAT: The code uses makesubimage, stiff, and dfits_theli programs. Without having those programs in set in your computer
            it wont work.

Required inputs (-r) image file for r-channel, (-g) image file for g-channel, (-b) image file for b-channel
                (-o) output image name, [[-reg region file (circles only)], -con contour file]

EXAMPLES: 
-         python3 create_subimages -r r_img.fits -g g_img.fits -b i_img.fits -o myimg -reg myreg.reg
          creates an image from the rgi images with maximum posibble resolution. If the region file contains single region
          The images are created only covering that region. The result will be myimg.tif and adjusted_myimg.tif with coordinates and region overlaid.
          If the region file has several region, image is created from the all area of the input images and regions will be overlaid. 

-         python3 create_subimages -r r_img.fits -g g_img.fits -b i_img.fits -o myimg -reg myreg.reg -bin 5
          does the same as in previous example. Image resolution is decreased by 5 by binning every 5 pixels.

-         python3 create_subimages -r r_img.fits -g g_img.fits -b i_img.fits -o myimg -reg myreg.reg -con ds9.con
          The given contour file is overlaid on the adjusted output image.

AUTHOR:
Melih Kara (karamel@astro.uni-bonn)


"""
)
parser.add_argument('-r', '--r_img', nargs=1,
                    required=True, help='input r-image name')
parser.add_argument('-g', '--g_img', nargs=1,
                    required=True, help='input g-image name')
parser.add_argument('-b', '--b_img', nargs=1,
                    required=True, help='input b-image name')
parser.add_argument('-o', '--output_file', nargs=1,
                    required=True, help='output file ')
parser.add_argument('-reg', '--region_file', nargs=1, default=['None'],required=False,
                    type=str, help='region file (OPTIONAL) \n ex: -reg ds9.reg')
parser.add_argument('-reg2', '--region_file2', nargs=1, default=['None'],required=False,
                    type=str, help='2nd region file (OPTIONAL) \n ex: -reg2 ds9_extra.reg')
parser.add_argument('-con', '--contour_file', nargs=1, required=False, default=['None'],
                    type=str, help='contour file\n ex: -con ds9.con')
parser.add_argument('-bin', '--binning', nargs=1, default=[1],
                    type=int, help='Choose a binning, default is 1. For large data resolution \
                                    can be decreased by binning the data (OPTIONAL) \n ex: -bin 5')
parser.add_argument('-xray', '--xray', nargs=1, default=[0],
                    type=int, help='If the image(s) are x-ray image set to 1 for proper images. (OPTIONAL)')

args = parser.parse_args()

# give meaningful variablenames to the command line
# arguments:
r_img = args.r_img[0]
g_img = args.g_img[0]
b_img = args.b_img[0]
output_fig = args.output_file[0]
regionfile = args.region_file[0]
reg2 = args.region_file2[0]
input_cont = args.contour_file[0]
binning = args.binning[0]
xray = args.xray[0]


# Initially consider to create an image from the all area of the fits. And overlay nothing
allimage = True ; regoverlay = False

def make_img(img1,img2,img3,outputname, overlay,bins=1):
    """
    Takes 3 fits images in i-r-g order, outputname. 
    CAVEAT: Runs *stiff* program on the terminal
    Creates one tif image and one png image.
    Optional binning can be specified to relax the process (less quality image)
    """
    if xray == 0:
        sp.call(["stiff", f"{img1}", f"{img2}", f"{img3}", "-MIN_LEVEL", "-0.2", "-MAX_LEVEL", "10.0", "-SKY_TYPE", "MANUAL", "-SKY_LEVEL", "0.0,0.0,0.0", \
                     "-BINNING", f"{bins}", "-MIN_TYPE", "MANUAL", "-MAX_TYPE", "MANUAL", "-OUTFILE_NAME", f"{outputname}.tif"]);
    else:
        sp.call(["stiff", f"{img1}", f"{img2}", f"{img3}","-c","xray.conf","-BINNING", f"{bins}", "-OUTFILE_NAME", f"{outputname}.tif"]);
    # call another python script to adjust the image
    if overlay:
        if reg2 == 'None':
            sp.call(["python3","adjust_image.py","-img",f"{outputname}.tif","-f",f"{img1}","-o", f"adjusted_{outputname}.tif",\
                     "-reg", f"{regionfile}", "-bin",f"{bins}"])
        else:
            sp.call(["python3","adjust_image.py","-img",f"{outputname}.tif","-f",f"{img1}","-o", f"adjusted_{outputname}.tif",\
                     "-reg", f"{regionfile}", "-bin",f"{bins}","-reg2", f"{reg2}" ])
    else:
        sp.call(["python3","adjust_image.py","-img", f"{outputname}.tif","-f",f"{img1}","-o", f"adjusted_{outputname}.tif", "-bin",f"{bins}"])
    print("IMAGE ADJUSTED")

# check if any region is given
if regionfile != 'None':
    try:
        regf = open(regionfile)
        reglines = regf.readlines()
        # Making sure that there are no comments or ds9 features in the beginning of the file
        ind = len([i for i in reglines[0:10] if 'circle' not in i])
        reglines = reglines[ind:]
        if len(reglines) == 1:           # If there is a single region 
            allimage = False             # do not create an image from the whole area
        regoverlay = True                # overlay the region in both cases
    except:
        print("Region file not found!")

if allimage == True:
    # call the script to create and overlay image
    make_img(b_img, r_img, g_img, output_fig, overlay=regoverlay, bins=binning)
    sys.exit(1)
else:
    singleregimage()

def singleregimage():
    '''
    Create an image around the given, single region
    Returns: output tif images, raw and adjusted. 
    '''
    # If we are not creating an image from the whole area
    # We have to determine the borders of the region file
    center_circ = [i.split('(')[1].split(',')[0:2] for i in reglines if i.split('(')[0]=='circle'][0]
    rad_circ = [i.split(')')[0].split(',')[-1] for i in reglines if i.split('(')[0]=='circle'][0]

    # If the radius is in arcsec convert to degrees
    try:
        rad_circ = float([rad_circ.split('"')[0]][0])/3600.
    except:
        print('The region file has to be in fk5 coordinates\n\
               The radii can either be in degrees or in arcsec with "')
        sys.exit(1)

    ra_deg = float(center_circ[0])
    dec_deg = float(center_circ[1])

    # Positions of the Lower Left corners in degrees
    RA_LL = ra_deg + rad_circ/np.cos(np.deg2rad(dec_deg))    # this is due to general image production convention i.e. RA decreases rightwards
    DEC_LL = dec_deg - rad_circ

    # Positions of the upper right corners in degrees
    RA_UR = ra_deg - rad_circ/np.cos(np.deg2rad(dec_deg))
    DEC_UR = dec_deg + rad_circ

    # We need to read the number of pixel in the image
    r_stream = os.popen(f"dfits_theli {r_img}  | grep NAXIS")                 # This assumes the number of pixels are stored in NAXIS column
    r_out = r_stream.read()                                                   # For an unclear reason, once the stream is read it loses the information
    r_n1 = r_out.split('=')[2].split('/')[0].strip()
    r_n2 = r_out.split('=')[3].split('/')[0].strip()
    g_stream = os.popen(f"dfits_theli {g_img}  | grep NAXIS")                
    g_out = g_stream.read()
    g_n1 = g_out.split('=')[2].split('/')[0].strip()
    g_n2 = g_out.split('=')[3].split('/')[0].strip()
    b_stream = os.popen(f"dfits_theli {b_img}  | grep NAXIS")                
    b_out = b_stream.read()
    b_n1 = b_out.split('=')[2].split('/')[0].strip()
    b_n2 = b_out.split('=')[3].split('/')[0].strip()

    # Define transformation functions
    def get_wcs_coor(filename, n1, n2):
        w = WCS(filename)
        min_RA, max_DEC = w.wcs_pix2world(n1, n2, 0)
        max_RA, min_DEC = w.wcs_pix2world(0, 0, 0)
        return float(min_RA), float(max_RA), float(min_DEC), float(max_DEC)

    def get_pix_coor(filename,RA_deg, DEC_deg):
        w = WCS(filename)
        RA_pix, DEC_pix = w.wcs_world2pix(RA_deg, DEC_deg, 0)
        return int(RA_pix), int(DEC_pix)

    # get the pixel coordinates, 
    # this ensures even if the orientation of different bands do not match
    # the pixel positions will correspond to correct object position
    rpix_ra_LL, rpix_dec_LL = get_pix_coor(r_img,RA_LL,DEC_LL)
    rpix_ra_UR, rpix_dec_UR = get_pix_coor(r_img,RA_UR,DEC_UR)
    rdx, rdy = abs(rpix_ra_UR-rpix_ra_LL) , abs(rpix_dec_UR-rpix_dec_LL)
    gpix_ra_LL, gpix_dec_LL = get_pix_coor(g_img,RA_LL,DEC_LL)
    gpix_ra_UR, gpix_dec_UR = get_pix_coor(g_img,RA_UR,DEC_UR)
    gdx, gdy = abs(gpix_ra_UR-gpix_ra_LL) , abs(gpix_dec_UR-gpix_dec_LL)
    bpix_ra_LL, bpix_dec_LL = get_pix_coor(b_img,RA_LL,DEC_LL)
    bpix_ra_UR, bpix_dec_UR = get_pix_coor(b_img,RA_UR,DEC_UR)
    bdx, bdy = abs(bpix_ra_UR-bpix_ra_LL) , abs(bpix_dec_UR-bpix_dec_LL)

    # Create a cut-out around the circle region
    # SYNTAX: makesubimage LowerLeftRAPixel LowerLeftDECPixel dx dy < fits_file >  OutputFitsFile
    #os.popen(f"makesubimage {rpix_ra_LL} {rpix_dec_LL} {rdx} {rdy} < {r_img} > cutout_{r_img}");
    #os.popen(f"makesubimage {gpix_ra_LL} {gpix_dec_LL} {gdx} {gdy} < {g_img} > cutout_{g_img}");
    #os.popen(f"makesubimage {bpix_ra_LL} {bpix_dec_LL} {bdx} {bdy} < {b_img} > cutout_{b_img}");

    fil = open("temp.sh","w")
    fil.write(f"makesubimage {rpix_ra_LL} {rpix_dec_LL} {rdx} {rdy} < {r_img} > cutout_{r_img}\n")
    fil.write(f"makesubimage {gpix_ra_LL} {gpix_dec_LL} {gdx} {gdy} < {g_img} > cutout_{g_img}\n")
    fil.write(f"makesubimage {bpix_ra_LL} {bpix_dec_LL} {bdx} {bdy} < {b_img} > cutout_{b_img}")
    fil.close()
    sp.run(["bash","temp.sh"])
    sp.run(["rm","temp.sh"])
    # call the make_img function to create image from the cutout fits.
    make_img(f'cutout_{b_img}', f'cutout_{r_img}', f'cutout_{g_img}', output_fig, overlay=regoverlay, bins=binning)

 



