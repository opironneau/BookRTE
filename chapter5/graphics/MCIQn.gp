# Monte Carlo with a variable refractive index (chapters 5 and 6)
set terminal pdfcairo size 12in,7.6in font 'Helvetica,13'
set output 'MCIQn.pdf'
set multiplot layout 2,2
  set xlabel 'z'; set ylabel 'T  [Celsius]'; set grid
  set key bottom left spacing 1.2 samplen 3; set xrange [0:1.2]
  set title 'temperature'
  set arrow 1 from 0.6,-32 to 0.6,20 nohead lc rgb 'gray50' dt 3
  set label 1 'Y' at 0.615,18 tc rgb 'gray30'
  plot 'cmp5_jump.txt'   u 1:2 w l lw 3 lc rgb '#d62728' t 'n = 1 | 0.9  (jump at Y)', \
       'cmp5_graded.txt' u 1:2 w l lw 3 lc rgb '#ff7f0e' t 'n = 1 {/Symbol \256} 0.9 (graded)', \
       'cmp5_flat.txt'   u 1:2 w l lw 3 lc rgb '#1f77b4' t 'n = 1', \
       'cmp5_sn.txt'     u 1:2 w l lw 2 dt 2 lc rgb 'black' t 'discrete ordinates, n = 1'
  unset arrow 1; unset label 1
  set title 'angular distribution on both sides of Y'
  set xlabel '{/Symbol m}'; set ylabel 'I/n^2'; set xrange [-1:1]; set autoscale y
  set key top left
  set arrow 2 from 0.4359,graph 0 to 0.4359,graph 1 nohead lc rgb 'gray40' dt 2
  set arrow 3 from -0.4359,graph 0 to -0.4359,graph 1 nohead lc rgb 'gray40' dt 2
  set label 2 '{/Symbol \261}{/Symbol m}_c' at 0.46,graph 0.93 tc rgb 'gray30'
  plot 'cmp5_ang_below.txt' u 1:2 w l lw 3 lc rgb '#d62728' t 'just below Y  (n=1)', \
       'cmp5_ang_above.txt' u 1:2 w l lw 3 lc rgb '#1f77b4' t 'just above Y  (n=0.9)'
  unset arrow 2; unset arrow 3; unset label 2
  set title 'I/n^2 (z,{/Symbol m})'
  set xlabel 'z'; set ylabel '{/Symbol m}'; set yrange [-1:1]; set xrange [0:1.2]
  set view map; set pm3d interpolate 2,2; unset key
  set arrow 1 from 0.6076,-1,9e9 to 0.6076,1,9e9 nohead lc rgb 'white' lw 2 front
  set arrow 2 from 0,0.4359,9e9 to 1.2,0.4359,9e9 nohead lc rgb 'white' lw 1 dt 2 front
  set arrow 3 from 0,-0.4359,9e9 to 1.2,-0.4359,9e9 nohead lc rgb 'white' lw 1 dt 2 front
  splot 'MCIQn_I.txt' using 1:2:3 with pm3d notitle
  unset arrow 1; unset arrow 2; unset arrow 3
  unset pm3d; set view; set autoscale
  set title 'polarization of the light leaving the atmosphere'
  set xlabel '{/Symbol m}'; set ylabel '-Q/I'; set xrange [0:1]; set grid; set key top right
  plot 'MCIQn_exit.txt' u 1:4 w lp pt 7 ps 0.6 lw 2 lc rgb '#d62728' t 'n = 1 | 0.9'
unset multiplot
