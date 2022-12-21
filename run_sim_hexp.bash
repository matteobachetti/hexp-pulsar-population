#!/bin/bash

xmldir=$SIXTE/share/sixte/instruments/hex-p

let_xml=${xmldir}/let/hexp_let_HEO_ff.xml
het_xml=${xmldir}/het/hexp_het_ff.xml

merged_sinput=$1
exposure=$2
RA=`fstatistic $merged_sinput+1 RA 1 | head -1 | awk '{print $8}'`
Dec=`fstatistic $merged_sinput+1 Dec 1 | head -1 | awk '{print $8}'`

root=`basename $merged_sinput .simput`

$SIXTE/bin/runsixt \
     XMLFile=${let_xml} \
     RA=$RA Dec=$Dec \
     Prefix=${root}_let_ \
     Simput=$merged_sinput \
     Exposure=$exposure

$SIXTE/bin/imgev \
     EvtFile=${root}_let_evt.fits \
     Image=${root}_let_img.fits \
     CoordinateSystem=0 Projection=TAN \
     CRVAL1=$RA CRVAL2=$Dec \
     NAXIS1=512 NAXIS2=512 CRPIX1=256.5 CRPIX2=256.5\
     CDELT1=-3.72422566e-04 CDELT2=3.72422566e-04 \
     CUNIT1=deg CUNIT2=deg \
     history=true clobber=yes

# run twice for HET!
for i in 1 2; do

    $SIXTE/bin/runsixt \
         XMLFile=${het_xml} \
	 RA=$RA Dec=$Dec \
         Prefix=${root}_het_${i}_ \
         Simput=$merged_sinput \
         Exposure=$exposure

    $SIXTE/bin/imgev \
         EvtFile=${root}_het_${i}_evt.fits \
         Image=${root}_het_${i}_img.fits \
         CoordinateSystem=0 Projection=TAN \
         NAXIS1=640 NAXIS2=640 CRPIX1=320.5 CRPIX2=320.5\
         CDELT1=-3.46066508e-04 CDELT2=3.46066508e-04 \
         CUNIT1=deg CUNIT2=deg \
	 CRVAL1=$RA CRVAL2=$Dec \
         history=true clobber=yes

done
