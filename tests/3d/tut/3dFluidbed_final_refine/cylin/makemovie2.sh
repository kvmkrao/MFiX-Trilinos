#!/bin/bash

textm="MFiX\n\nNETL"
textt="With Trilinos\n\nSandia\n\n"
posm=+50+300
post=+0+300
delay=30
path1=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final_refine/cylin/MFIX/reduceddim_final/animation/
path2=/home/vkotteda/DOE/mfix_mpi3d/tut/3dFluidbed_final_refine/cylin/TRI/reduceddim_final/animation/
#text1="~50\\% \n faster "
text1=">50\\%  faster "
pos=+0+500


mkdir -p movie2
for i in 00 01 02 03 04 05 06 07 08 09 ; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix3d_cylepg.00$i.png  movie2/mfix3ddem_epg.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri3dcyl_epg.00$i.png   movie2/auxilary$i.png
 convert  -pointsize 32 -annotate $pos  "$text1" movie2/auxilary$i.png          movie2/tri3d_dem_epg.00$i.png
#  convert  -pointsize 32 -annotate $post  "$text1" movie2/tri2d_dem_epg.00$i.png  movie1/tri2d_dem_epg.00$i.png
 convert movie2/mfix3ddem_epg.00$i.png movie2/tri3d_dem_epg.00$i.png +append    movie2/mfixtri3d_dem_epg.00$i.png
 convert movie2/mfixtri3d_dem_epg.00$i.png    ~/Documents/names8.png -append    movie2/final.00$i.png
convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" movie2/final.00$i.png movie2/finaln.00$i.png
done

for i in $(seq 10 1 60) ; do
 echo $i
 convert  -pointsize 32 -annotate $posm "$textm" $path1/mfix3d_cylepg.00$i.png     movie2/mfix3ddem_epg.00$i.png
 convert  -pointsize 32 -annotate $post "$textt" $path2/tri3dcyl_epg.00$i.png      movie2/auxilary$i.png
 convert  -pointsize 32 -annotate $pos  "$text1" movie2/auxilary$i.png             movie2/tri3d_dem_epg.00$i.png
#  convert  -pointsize 32 -annotate $post  "$text1" movie2/tri2d_dem_epg.00$i.png  movie1/tri2d_dem_epg.00$i.png
 convert movie2/mfix3ddem_epg.00$i.png movie2/tri3d_dem_epg.00$i.png +append       movie2/mfixtri3d_dem_epg.00$i.png
 convert movie2/mfixtri3d_dem_epg.00$i.png    ~/Documents/names8.png -append       movie2/final.00$i.png
 convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" movie2/final.00$i.png movie2/finaln.00$i.png 
done

(cd movie2; convert -delay $delay -quality 100 -normalize finaln.???[1,3,5,7,9].png movie2/final.0000.png movie2/final.0000.png finaln.???[1,3,5,7,9].png mfixtri3dcyltfm_epg.gif; rm mfix3ddem_epg* tri3d_dem_epg*  mfixtri3d_dem_epg* auxilary* )
(cd movie2; convert -delay $delay -quality 100 -normalize finaln.???[1,3,5,7,9].png  movie2/final.0000.png movie2/final.0000.png finaln.???[1,3,5,7,9].png mfixtri3dcyltfm_epg.mp4)
