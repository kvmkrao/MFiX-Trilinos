#!/bin/bash

textm="MFiX\n\n"
textt="MFiX-Trilinos\n\n"
posm=+50+300
post=+0+300
delay=50
path1=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/MFIX/animate
path2=/home/vkotteda/DOE/mfix_mpi2d/tut/FluidBed_DES/TRI/animation
text1="~50\\% faster "
pos=+0+500

dir=$(pwd)/movie_clean
mkdir -p $dir
for i in 00 02 04 06 08 10 12 14 16 18 20; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix_dem2d.00$i.png     $dir/mfix2ddem_epg.$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri2d_dem_epg.00$i.png  $dir/auxilary$i.png
 convert $dir/mfix2ddem_epg.$i.png $dir/auxilary$i.png +append    $dir/final.$i.png
done

(cd $dir; convert -delay $delay -quality 100 -normalize  final.00.png final.??.png mfixtri2ddem_epg.gif; rm mfix2ddem_epg* tri2d_dem_epg* aux*)
#(cd $dir; convert -delay $delay -quality 100 -normalize  finaln.????.png final.0000.png finaln.????.png final.0000.png finaln.????.png mfixtri2ddem_epg.mp4)
