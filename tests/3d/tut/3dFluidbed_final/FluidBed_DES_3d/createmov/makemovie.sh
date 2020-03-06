#!/bin/bash

textm="MFiX\n\nNETL"
textt="With Trilinos\n\nSandia\n\n"
posm=+50+300
post=+0+300
delay=50
path1=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final/FluidBed_DES_3d/MFIX
path2=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final/FluidBed_DES_3d/TRI
#text1="~50\\% \n faster "
text1=">50\\%  faster "
pos=+0+500


mkdir -p movie1
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 ; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix_dem3d.00$i.png     movie1/mfix.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri3d_epg.00$i.png      movie1/auxilary$i.png
 convert  -pointsize 32 -annotate $pos  "$text1" movie1/auxilary$i.png          movie1/tri.00$i.png
# convert  -pointsize 32 -annotate $post  "$text1" movie1/tri2d_dem_epg.00$i.png  movie1/tri2d_dem_epg.00$i.png
  convert movie1/mfix.00$i.png movie1/tri.00$i.png   +append                    movie1/auxilary1.00$i.png
  convert movie1/auxilary1.00$i.png  ~/Documents/names8.png   -append      movie1/final.00$i.png 
  convert -pointsize 32 -annotate +250+50 "Fluidized Bed (DEM, gas-solid)" movie1/final.00$i.png movie1/finaln.00$i.png
done

(cd movie1; convert -delay $delay -quality 100 -normalize finaln.????.png mfixtri3ddem_epg.gif) #; rm mfix3ddem_epg* tri3d_dem_epg* mfix3ddem_epg*   auxilary* auxilary1*)
