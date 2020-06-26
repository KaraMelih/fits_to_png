#!/bin/sh

python3 create_subimages.py -r test/test_r_imgs.fits -g test/test_g_imgs.fits -b test/test_i_imgs.fits \
                            -o test/output/test1

python3 create_subimages.py -r test/test_r_imgs.fits -g test/test_g_imgs.fits -b test/test_i_imgs.fits \
                            -o test/output/test2 -con test/3Gsm_6lvl9sm.con

python3 create_subimages.py -r test/test_r_imgs.fits -g test/test_g_imgs.fits -b test/test_i_imgs.fits \
                            -o test/output/test3 -reg test/test.reg

python3 create_subimages.py -r test/testXray_b3.fits -g test/testXray_b2.fits -b test/testXray_b8.fits \
                            -o test/output/test4_xray -xray 1

python3 create_subimages.py -r test/test_r_imgs.fits -g test/test_g_imgs.fits -b test/test_i_imgs.fits \
                            -o test/output/test5 -reg test/test2.reg

