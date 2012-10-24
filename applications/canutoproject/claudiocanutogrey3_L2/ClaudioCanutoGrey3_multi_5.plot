reset

set output "ClaudioCanutoGrey3_multi_5.png" 
set size 1,1

set term png transparent size 4000,2453 crop


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


set origin
set multiplot

set size 0.5,1
set origin 0,0


splot 'image_0.7_0.1_standard_5.dat' u 1:2:3 w pm3d notitle

set size 0.5,1
set origin 0.5,0


splot 'image_0.7_0.1_refsol.dat' u 1:2:3 w pm3d notitle,\
      'image_coeff_0.7_0.1_standard_5.dat' u 1:2:3 w p ps 2 lc rgb 'red' notitle


set size 0.25,0.5
set origin 0.41,0.38

set xrange[350:1150]
set yrange[650:1250]

splot 'image_0.7_0.1_refsol.dat' u 1:2:3 w pm3d notitle,\
      'image_coeff_0.7_0.1_standard_5.dat' u 1:2:3 w p ps 2 lc rgb 'red' notitle


set size 0.25,0.5
set origin 0.35,0.11

set xrange[350:1150]
set yrange[650:1250]

splot 'image_0.7_0.1_standard_5.dat' u 1:2:3 w pm3d notitle

unset multiplot




