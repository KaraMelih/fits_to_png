The script creates single color or three color images from the given fits files.
If ds9 region file with the circle region(s) is given.
If there is a single region, it creates the image around this given region and also draws the cirle around it. 
If there are couple of regions in the file, it displays all the region on the whole image. 
It can also take ds9 contour files and overlays on the images.

Default is for optical images, but one can also create x-ray images using x-ray configuration file.

CAVEATS: It uses, makesubimage, stiff and swarp programs. :(

- **It is a work in progress so there might be some issues**
