#!/bin/bash

textm="MFiX\n\n"
textt="MFiX-Trilinos\n\n"
posm=+150+200
post=+480+200
delay=30
path1=/home/vkotteda/DOE/mfix_mpi2d/tut/fluidBed1_refine//MFIX/animation1
path2=/home/vkotteda/DOE/mfix_mpi2d/tut/fluidBed1_refine/TRI
#text1="~50\\% \n faster "
text1="~50\\%  faster "
pos=+450+500


dir=$(pwd)/movie_plain
mkdir -p $dir

for i in 00 01 02 03 04 05 06 07 08 09 ; do
 echo $i

  convert   $path1/epg_mfix2d.00$i.png   $path2/epg_tri2d.00$i.png +append  $dir/auxilary$i.png 
#  convert   $dir/auxilary$i.png ~/Documents/names8_resize.png                          -append  $dir/final1.00$i.png
  convert  -pointsize 32 -annotate $posm "$textm" $dir/auxilary$i.png   $dir/final2.00$i.png
  convert  -pointsize 32 -annotate $post "$textt" $dir/final2.00$i.png  $dir/final3.00$i.png
#  convert  -pointsize 32 -annotate $pos  "$text1" $dir/final3.00$i.png  $dir/final.00$i.png
#  convert  -pointsize 32 -annotate +250+30 "Fluidized Bed (TFM, gas-soild)" $dir/final.00$i.png  $dir/finaln.00$i.png
done

for i in $(seq 10 1 99) ; do
 echo $i
    convert $path1/epg_mfix2d.00$i.png $path2/epg_tri2d.00$i.png +append  $dir/auxilary$i.png
#  convert $dir/auxilary$i.png ~/Documents/names8_resize.png                          -append  $dir/final1.00$i.png
  convert  -pointsize 32 -annotate $posm "$textm" $dir/auxilary$i.png                        $dir/final2.00$i.png
  convert  -pointsize 32 -annotate $post "$textt" $dir/final2.00$i.png                       $dir/final3.00$i.png
#  convert  -pointsize 32 -annotate $pos  "$text1" $dir/final3.00$i.png                       $dir/final.00$i.png
#  convert  -pointsize 32 -annotate +250+30 "Fluidized Bed (TFM, gas-soild)" $dir/final.00$i.png  $dir/finaln.00$i.png
done

(cd $dir; convert -delay $delay -quality 100 -normalize final3.0??[1,3,5,7,9].png mfix2dtfm_epg.gif);
(cd $dir; convert -delay $delay -quality 100 -normalize final3.0??[1,3,5,7,9].png mfix2dtfm_epg.mp4; rm final1* final2* auxilary* )
