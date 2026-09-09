# gnuplot milneMC_I.gp
# col=3: diffuse intensity, the direct beam being drawn apart as a Dirac;
# col=4: the beam is instead smeared over the mu bin which contains mu_s.
col = 3
mus = 0.57735
set terminal pngcairo size 1500,760 enhanced font 'Helvetica,13'
set output 'milneMC_I.png'
set palette defined (0 '#2b1b6e', 0.35 '#1f8ac0', 0.6 '#7fd06a', 0.8 '#f6d746', 1 '#b3202a')
set xlabel 'z'; set ylabel '{/Symbol m}'
set xrange [0:1]; set yrange [-1:1]
set multiplot layout 1,2 title "Milne problem by Monte Carlo: diffuse intensity I(z,{/Symbol m}),   Z=1, {/Symbol k}=1, {/Symbol s}_s=1, {/Symbol m}_s=0.57735,  2000000 photons"
  set pm3d depthorder
  set hidden3d; set style fill transparent solid 0.9
  set contour base; set cntrparam levels 12; unset clabel
  set zlabel 'I' rotate by 90
  set view 55,340,1,1.15
  set ticslevel 0.05
  unset colorbox; unset key
  set zrange [0:1.05]
  set label 11 'Dirac {/Symbol d}({/Symbol m}-{/Symbol m}_s), of weight e^{-{/Symbol k}z/{/Symbol m}_s}' at 0.05*1,mus,1.0 tc 'red' front
  set label 12 'crest = {/Symbol s}_sJ_0(z)/{/Symbol k} : the grazing limit I(z,0^{/*0.8 +}_{/*0.8 -})' at 0.45*1,0,0.80 tc 'black' front
  splot 'milneMC_I.txt' using 1:2:col with pm3d notitle, \
        'milneMC_I.txt' using 1:2:col with lines lc 'black' lw 0.2 nogrid notitle, \
        'milneMC_J0.txt' using 1:(0):5 with lines lw 3 lc 'black' nogrid nocontours notitle, \
        'milneMC_J0.txt' every 2 using 1:(mus):6 with impulses lw 1.5 lc 'red' nogrid nocontours notitle
  set colorbox; set view map; set pm3d interpolate 2,2
  unset label 11; unset label 12; unset zrange
  set arrow 1 from 0,mus,9e9 to 1,mus,9e9 nohead lc 'red' lw 2 dt 2 front
  set label 1 '  {/Symbol d}({/Symbol m}-{/Symbol m}_s), the incoming beam' at 0,mus,9e9 front tc 'red' offset 0,0.7
  set arrow 2 from 0,0,9e9 to 1,0,9e9 nohead lc 'black' lw 1 dt 3 front
  set label 2 '  {/Symbol m}=0: I = {/Symbol s}_sJ_0/{/Symbol k}' at 0.55*1,0,9e9 front offset 0,0.7
  set title 'seen from above, with the level lines'
  splot 'milneMC_I.txt' using 1:2:col with pm3d notitle, \
        'milneMC_I.txt' using 1:2:col with lines lc 'black' lw 0.6 nosurface notitle
unset multiplot
