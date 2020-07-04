The script creates single color or three color images from the given fits files.
If ds9 region file with the circle region(s) is given.
If there is a single region, it creates the image around this given region and also draws the cirle around it. 
If there are couple of regions in the file, it displays all the region on the whole image. 
It can also take ds9 contour files and overlays on the images.

Default is for optical images, but one can also create x-ray images using x-ray configuration file.

CAVEATS: It uses, makesubimage, stiff and swarp programs. :(

- **It is a work in progress so there might be some issues**
- currently cannot make single color image. 3 channels are required


# x-ray_color.py

**Can be used alone, syntax:** <br>

python3 x-ray_color.py -f fitsfile.fits -o outputname ((-cmap inferno) (-sm 3) (-cr True) (-br 1.7)) <br>

call ```python3 x-ray_color.py -h``` to display help. <br>

Uses <font color='red'> Stiff </font> Program from Astromatic ! <br>

Example output; grey scale tif image and colored image with the given parameters.
<img src="test.tif">
<img src="test.png">