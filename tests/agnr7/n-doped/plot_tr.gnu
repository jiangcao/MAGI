

#set term pngcairo size 500,400 enhanced font 'Arial, 14'
#set output 'fig_ldos_off.png'
unset key
load '~/gnuplot/gnuplot-palettes/parula.pal'
# set yrange [-4:0]
#set xrange [0:140]
#set cbrange [0:30]

set xlabel 'x (nm)'
set ylabel 'E (eV)'

id='00'
p \
 'bse_trL00'.id.'.dat' w l,\
 'bse_trR00'.id.'.dat' w l

pause -1

