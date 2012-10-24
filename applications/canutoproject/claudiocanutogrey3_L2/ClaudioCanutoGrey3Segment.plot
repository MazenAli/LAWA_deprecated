reset
set term png transparent size 601,801 crop
set output "ClaudioCanutoGrey3Segment_30.png" 

set view 0,90
set palette gray

set xrange[350:1150]
set yrange[650:1250]


unset border
unset colorbox
unset xlabel
unset xtics
unset ylabel
unset ytics
unset zlabel
unset ztics
unset clabel
unset cbtics
splot 'image_0.7_0.1_standard_30.dat' u 1:2:3 w pm3d notitle
#      'image_coeff_0.7_0.1_standard_30.dat' u 1:2:3 w p ps 2 lc rgb 'red' notitle


set output "ClaudioCanutoGrey3Segment_30_approx.png" 

set view 0,90
set palette gray

set xrange[350:1150]
set yrange[650:1250]


unset border
unset colorbox
unset xlabel
unset xtics
unset ylabel
unset ytics
unset zlabel
unset ztics
unset clabel
unset cbtics
splot 'image_0.7_0.1_standard_30.dat' u 1:2:4 w pm3d notitle
#      'image_coeff_0.7_0.1_standard_30.dat' u 1:2:3 w p ps 2 lc rgb 'red' notitle

