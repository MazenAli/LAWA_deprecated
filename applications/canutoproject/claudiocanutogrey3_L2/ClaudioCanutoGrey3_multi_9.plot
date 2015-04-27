reset

set output "ClaudioCanutoGrey3_multi_9.png" 
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
set origin 0.025,0


splot 'image_L2_0.7_0.1_standard_sparsetree__9.dat' u 1:2:3 w pm3d notitle

set size 0.5,1
set origin 0.475,0


splot 'image_L2_0.7_0.1_standard_sparsetree__refsol.dat' u 1:2:3 w pm3d notitle,\
      'scattercoeff_L2_0.7_0.1_standard_sparsetree__9.dat' u 1:2:3 w p ps 1 pt 3 lc rgb 'red' notitle


set size 0.25,0.5
set origin 0.375,0.39

set xrange[350:1150]
set yrange[650:1250]

splot 'image_L2_0.7_0.1_standard_sparsetree__refsol.dat' u 1:2:3 w pm3d notitle,\
      'scattercoeff_L2_0.7_0.1_standard_sparsetree__9.dat' u 1:2:3 w p ps 1 pt 3 lc rgb 'red' notitle


set size 0.25,0.5
set origin 0.375,0.11

set xrange[350:1150]
set yrange[650:1250]

splot 'image_L2_0.7_0.1_standard_sparsetree__9.dat' u 1:2:3 w pm3d notitle

unset multiplot