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


mkdir -p movie
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 ; do
 echo $i
  convert   $path1/mfix_dem3d.00$i.png $path2/tri3d_epg.00$i.png    +append      movie/auxilary1.00$i.png
done

(cd movie; convert -delay $delay -quality 100 -normalize auxilary1.????.png mfixtri3ddem_epg.gif; rm auxilary1.* ) 
