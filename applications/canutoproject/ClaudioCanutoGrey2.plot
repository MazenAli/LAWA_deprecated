reset
set term png transparent size 326,400 crop
set output "ClaudioCanutoGrey2_31.png" 

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
splot 'image_poisson_0.7_0.001_1e-06_standard_sparsetree_31.dat' u 1:2:3 w pm3d notitle
#      'scattercoeff_image_poisson_0.7_0.001_standard_sparsetree_30.dat' u 1:2:3 w p ps 0.1 lc rgb 'red' notitle


set output "ClaudioCanutoGrey2_approx_31.png" 

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
splot 'image_poisson_0.7_0.001_1e-06_standard_sparsetree_31.dat' u 1:2:4 w pm3d notitle

