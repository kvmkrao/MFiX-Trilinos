#!/bin/bash

textm="MFiX\n\n "
textt="MFiX-Trilinos \n\n"
posm=+50+300
post=+0+300
delay=30
path1=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final_refine/cart/MFIX/animation
path2=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final_refine/cart/TRI/animation
#text1="~50\\% \n faster "
text1=">50\\%  faster "
pos=+0+500

dir=$(pwd)/movie_plain

mkdir -p $dir
for i in 00 01 02 03 04 05 06 07 08 09 ; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix3d_epg.00$i.png     $dir/mfix3ddem_epg.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri3d_epg.00$i.png      $dir/auxilary$i.png
# convert  -pointsize 32 -annotate $pos  "$text1" $dir/auxilary$i.png          $dir/tri3d_dem_epg.00$i.png
#  convert  -pointsize 32 -annotate $post  "$text1" $dir/tri2d_dem_epg.00$i.png  movie1/tri2d_dem_epg.00$i.png
 convert $dir/mfix3ddem_epg.00$i.png $dir/auxilary$i.png +append    $dir/final.00$i.png
# convert $dir/mfixtri3d_dem_epg.00$i.png    ~/Documents/names8.png -append    $dir/final.00$i.png
#  convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" $dir/final.00$i.png      $dir/finaln.00$i.png
done


for i in $(seq 10 1 99) ; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix3d_epg.00$i.png     $dir/mfix3ddem_epg.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri3d_epg.00$i.png      $dir/auxilary$i.png
# convert  -pointsize 32 -annotate $pos  "$text1" $dir/auxilary$i.png          $dir/tri3d_dem_epg.00$i.png
 convert $dir/mfix3ddem_epg.00$i.png $dir/auxilary$i.png +append  $dir/final.00$i.png
# convert $dir/mfixtri3d_dem_epg.00$i.png    ~/Documents/names8.png -append    $dir/final.00$i.png
#  convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" $dir/final.00$i.png      $dir/finaln.00$i.png
done

(cd $dir; convert -delay $delay -quality 100 -normalize final.???[1,3,5,7,9].png final.0000.png final.0000.png final.???[1,3,5,7,9].png mfixtri3dtfm_epg.gif; rm mfix3ddem_epg* tri3d_dem_epg* mfixtri3d_dem_epg* auxilary* )

#(cd $dir; convert -delay $delay -quality 100 -normalize finaln.???[1,3,5,7,9].png final.0000.png final.0000.png finaln.???[1,3,5,7,9].png mfixtri3dtfm_epg.gif)
