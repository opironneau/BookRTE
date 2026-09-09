# gnuplot MCIQ4.gp
set terminal pngcairo size 1350,900 enhanced font 'Helvetica,13'
set output 'MCIQ4.png'
set palette defined (0 '#2b1b6e', 0.35 '#1f8ac0', 0.6 '#7fd06a', 0.8 '#f6d746', 1 '#b3202a')
mus=0.5
set multiplot layout 2,2 title "Monte Carlo for (4.9)-(4.11):  {/Symbol b}=0.5, {/Symbol m}_s=0.5, q_0=-0.3, c'_e=2, c_s=2e-06;   surfaces at {/Symbol n}=0.678537 ({/Symbol l}=4.42128 {/Symbol m}m, {/Symbol k}_{/Symbol n}=0.5) with 4000000 photons"
  set xlabel 'z'; set ylabel '{/Symbol m}'; set xrange [0:0.696]; set yrange [-1:1]
  set pm3d depthorder; set hidden3d; set style fill transparent solid 0.9
  set ticslevel 0.05; set view 55,340,1,1.15; unset key; unset colorbox
  set title 'diffuse intensity I(z,{/Symbol m}); red: the Dirac of the sunlight'
  splot 'MCIQ4_I.txt' using 1:2:3 with pm3d notitle, \
        'MCIQ4_I.txt' using 1:2:3 with lines lc 'black' lw 0.2 nogrid nocontours notitle, \
        'MCIQ4_Jnu.txt' every 2 using 1:(-mus):6 with impulses lw 1.5 lc 'red' \
             nogrid nocontours notitle
  set title 'polarization Q(z,{/Symbol m})'
  splot 'MCIQ4_I.txt' using 1:2:4 with pm3d notitle, \
        'MCIQ4_I.txt' using 1:2:4 with lines lc 'black' lw 0.2 nogrid nocontours notitle
  unset pm3d; unset hidden3d; set view map; set grid; set key top right
  set title 'temperature, solution of (4.10)'
  set xlabel 'z'; set ylabel 'T  [Celsius]'; set autoscale x; set autoscale y
  plot 'MCIQ4_T.txt' using 1:2 with lines lw 2 lc 'blue' title 'T(z)'
  set title 'degree of polarization of the light leaving the atmosphere'
  set xlabel '{/Symbol m}'; set ylabel '-Q(Z,{/Symbol m})/I(Z,{/Symbol m})'
  set xrange [0:1]
  plot 'MCIQ4_exit.txt' using 1:4 with linespoints pt 7 ps 0.6 lw 2 lc 'red' \
       title 'at {/Symbol n}=0.678537'
unset multiplot
