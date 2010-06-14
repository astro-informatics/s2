#!/bin/sh

rm -f test/gcmb.*

bin/s2_gcmb -inp test/wmap_lcdm_pl_model_yr1_v1.txt \
  -scale_cl .true. -out test/gcmb.fits -nside 128 -lmax 256 -seed 1

map2gif -inp test/gcmb.fits -out test/gcmb.gif -bar .true.

