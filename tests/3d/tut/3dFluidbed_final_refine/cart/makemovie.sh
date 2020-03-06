
#!/bin/bash

textm="MFIX\n\nNETL"
textt="With Trilinos\n\nSandia\n\n"
posm=+50+300
post=+450+300
delay=20
path1=/home/vkotteda/DOE/mfix_mpi2d/tut/fluidBed1_refine/increased_width/MFIX/animation1
path2=/home/vkotteda/DOE/mfix_mpi2d/tut/fluidBed1_refine/increased_width/TRI
#text1="~50\\% \n faster "
text1="~50\\%  faster "
pos=+450+500


mkdir -p movie1

for i in 00 01 02 03 04 05 06 07 08 09 ; do
 echo $i

  convert blank.png $path1/epg_mfix2d.00$i.png blank.png  $path2/epg_tri2d.00$i.png +append  movie1/auxilary$i.png
  convert movie1/auxilary$i.png names3.png                          -append  movie1/final1.00$i.png
  convert  -pointsize 32 -annotate $posm "$textm" movie1/final1.00$i.png  movie1/final2.00$i.png
  convert  -pointsize 32 -annotate $post "$textt" movie1/final2.00$i.png  movie1/final3.00$i.png
  convert  -pointsize 32 -annotate $pos  "$text1" movie1/final3.00$i.png  movie1/final.00$i.png
  convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" movie2/final.00$i.png      movie2/finaln.00$i.png
done

for i in $(seq 10 1 99) ; do
 echo $i
    convert blank.png $path1/epg_mfix2d.00$i.png blank.png  $path2/epg_tri2d.00$i.png +append  movie1/auxilary$i.png
  convert movie1/auxilary$i.png names3.png                          -append  movie1/final1.00$i.png
  convert  -pointsize 32 -annotate $posm "$textm" movie1/final1.00$i.png  movie1/final2.00$i.png
  convert  -pointsize 32 -annotate $post "$textt" movie1/final2.00$i.png  movie1/final3.00$i.png
  convert  -pointsize 32 -annotate $pos  "$text1" movie1/final3.00$i.png  movie1/final.00$i.png
  convert -pointsize 26 -annotate +350+50 "Fluidized Bed (TFM, gas-solid)" movie2/final.00$i.png      movie2/finaln.00$i.png
done

(cd movie1; convert -delay $delay -quality 100 -normalize finaln.0??[1,3,5,7,9].png mfix3ddem_epg.gif; rm final1* final2* final3* )
(cd movie1; convert -delay $delay -quality 100 -normalize finaln.0??[1,3,5,7,9].png mfix3ddem_epg.mp4) 

