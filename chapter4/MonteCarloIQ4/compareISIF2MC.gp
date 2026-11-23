# Monte Carlo and ISIF, with and without CO2, all four on the same 120 point grid
set terminal pdfcairo size 9.2in,3.6in font 'Helvetica,13'
set output 'compareISIF2MC.pdf'
set multiplot layout 1,2
  set xlabel 'altitude z'; set ylabel 'T  [Celsius]'
  set grid; set key bottom left spacing 1.2 samplen 3
  set xrange [0:1.2]; set yrange [-45:5]
  set title 'temperature'
  plot 'cmp_mc_noco2.txt'   u 1:2 w l lw 3 lc rgb '#1f77b4' t 'Monte Carlo', \
       'cmp_mc_co2.txt'     u 1:2 w l lw 3 lc rgb '#d62728' t 'Monte Carlo, CO_2', \
       'cmp_isif_noco2.txt' u 1:2 w l lw 2 dt 2 lc rgb '#0b3d62' t 'ISIF', \
       'cmp_isif_co2.txt'   u 1:2 w l lw 2 dt 2 lc rgb '#7f1416' t 'ISIF, CO_2'
  set title 'effect of CO_2:  T(1.2 {/Symbol k}_{/Symbol n}) - T({/Symbol k}_{/Symbol n})'
  set ylabel '{/Symbol D}T  [Celsius]'; set autoscale y
  set key top right
  plot 0 lc rgb 'black' lw 1 notitle, \
       'cmp_dT_mc.txt'   u 1:2 w l lw 3 lc rgb '#1f77b4' t 'Monte Carlo', \
       'cmp_dT_isif.txt' u 1:2 w l lw 2 dt 2 lc rgb '#7f1416' t 'ISIF'
unset multiplot
