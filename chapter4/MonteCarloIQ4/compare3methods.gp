# the three methods on the same 119/120 point grid, CO2 case
set terminal pdfcairo size 9.2in,3.6in font 'Helvetica,13'
set output 'compare3methods.pdf'
set multiplot layout 1,2
  set xlabel 'altitude z'; set ylabel 'T  [Celsius]'
  set grid; set key bottom left spacing 1.2 samplen 3
  set xrange [0:1.2]; set yrange [-45:5]
  set title 'temperature, 554 frequencies, CO_2 case'
  plot 'cmp3_mc.txt'   u 1:2 w l lw 5 lc rgb '#2ca02c' t 'Monte Carlo', \
       'cmp3_sn.txt'   u 1:2 w l lw 2 lc rgb '#1f77b4' t 'discrete ordinates', \
       'cmp3_isif.txt' u 1:2 w l lw 2 dt 2 lc rgb '#d62728' t 'ISIF'
  set title 'difference with the discrete ordinates'
  set ylabel '{/Symbol D}T  [Celsius]'; set yrange [-0.8:0.8]
  set key top right
  plot 0 lc rgb 'black' lw 1 notitle, \
       'cmp3_dmc.txt'   u 1:2 w l lw 2 lc rgb '#2ca02c' t 'Monte Carlo', \
       'cmp3_disif.txt' u 1:2 w l lw 2 dt 2 lc rgb '#d62728' t 'ISIF'
unset multiplot
