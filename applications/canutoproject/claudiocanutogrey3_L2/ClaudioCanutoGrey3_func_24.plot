reset
set term png transparent size 2000,2453 crop

set output "ClaudioCanutoGrey3_func_24.png" 
set view 22,68
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
splot 'image_L2_0.7_0.1_standard_sparsetree__refsol.dat' u 1:2:3 w pm3d notitle


reset