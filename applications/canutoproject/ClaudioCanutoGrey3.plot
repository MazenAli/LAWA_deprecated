reset
set term png transparent size 2000,2453 crop
set output "ClaudioCanutoGrey3_7.png" 

set view 0,90
set palette gray


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
splot 'image_0.7_0.1_standard_7.dat' u 1:2:3 w pm3d notitle
#      'image_coeff_0.7_0.1_standard_30.dat' u 1:2:3 w p ps 2 lc rgb 'red' notitle

