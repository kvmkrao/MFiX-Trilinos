#!/bin/bash

textm="MFiX\n\nNETL"
textt="With Trilinos\n\nSandia\n\n"
posm=+50+300
post=+0+300
delay=20
path1=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/MFIX/animate
path2=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/TRI/animation
text1="~50\\% faster "
pos=+0+500

mkdir -p movie
for i in 00 02 04 06 08 10 12 14 16 18 20; do
 echo $i

# convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix_dem2d.00$i.png     movie1/mfix2ddem_epg.00$i.png
# convert  -pointsize 32 -annotate $post "$textt" $path2/tri2d_dem_epg.00$i.png  movie1/auxilary$i.png
# convert  -pointsize 32 -annotate $pos  "$text1" movie1/auxilary$i.png          movie1/tri2d_dem_epg.00$i.png
 convert    $path1/mfix_dem2d.00$i.png $path2/tri2d_dem_epg.00$i.png    +append  movie/auxilary1.00$i.png
# convert movie1/mfixtri2d_dem_epg.00$i.png names3.png -append  movie1/final.00$i.png
# convert movie1/auxilary1.00$i.png ~/Documents/names8.png -append  movie1/final.00$i.png
done

(cd movie; convert -delay $delay -quality 100 -normalize auxilary1.????.png mfixtri2ddem_epg.gif; rm auxilary1* )
