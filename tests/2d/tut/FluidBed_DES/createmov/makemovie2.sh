#!/bin/bash

textm="MFiX\n\nNETL"
textt="With Trilinos\n\nSandia\n\n"
posm=+50+300
post=+0+300
delay=50
path1=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/MFIX/animate
path2=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/TRI/animation
text1="~50\\% faster "
pos=+0+500

mkdir -p movie2
for i in 00 02 04 06 08 10 12 14 16 18 20; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix_dem2d.00$i.png     movie2/mfix2ddem_epg.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri2d_dem_epg.00$i.png  movie2/auxilary$i.png
 convert  -pointsize 32 -annotate $pos  "$text1" movie2/auxilary$i.png          movie2/tri2d_dem_epg.00$i.png
 convert movie2/mfix2ddem_epg.00$i.png movie2/tri2d_dem_epg.00$i.png +append    movie2/auxilary1.00$i.png
# convert movie2/mfixtri2d_dem_epg.00$i.png names3.png -append  movie2/final.00$i.png
 convert movie2/auxilary1.00$i.png ~/Documents/names8.png -append  movie2/final.00$i.png
 convert -pointsize 26 -annotate +350+50 "Fluidized bed (DEM, gas-soild)" movie2/final.00$i.png  movie2/finaln.00$i.png
done

(cd movie2; convert -delay $delay -quality 100 -normalize finaln.????.png final.0000.png finaln.????.png final.0000.png finaln.????.png mfixtri2ddem_epg.gif; rm mfix2ddem_epg* tri2d_dem_epg* aux*)
(cd movie2; convert -delay $delay -quality 100 -normalize finaln.????.png final.0000.png finaln.????.png final.0000.png finaln.????.png mfixtri2ddem_epg.mp4)
