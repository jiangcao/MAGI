

#set term pngcairo size 500,400 enhanced font 'Arial, 14'
#set output 'fig_ldos_off.png'
set pm3d map
unset key
load '~/gnuplot/gnuplot-palettes/parula.pal'
set yrange [-4:0]
#set xrange [0:140]
set cbrange [0:30]

set xlabel 'x (nm)'
set ylabel 'E (eV)'


id='00'
sp \
 'gw_ldos00'.id.'.dat'

pause -1


id='01'
sp \
 'gw_ldos00'.id.'.dat'

pause -1

