set pm3d map                                                                                                      
#set yrange [-5.5:1]
#set xrange [-1:23]
set cbrange [0:30]
unset key 
load '/home/jiacao/gnuplot/gnuplot-palettes/moreland.pal'

set multiplot layout 1,3

dir='./'

iter='0000'
sp  dir.'gw_ldos'.iter.'.dat'

iter='0001'
sp  dir.'gw_ldos'.iter.'.dat'

iter='0000'
sp  dir.'bse_ldos'.iter.'.dat'

unset multiplot

pause -1
