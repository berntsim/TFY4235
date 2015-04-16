#reset
#set term gif animate
#set output "testAnimation.gif"
#n = 24  
#dt = 2*pi/n
#set xrange [0:1]
#set yrange [0:1]
#do for [i=0:n] {
#        splot
#}
set xlabel "x-points"
set ylabel "y-points"

set view 60,45
set pm3d
splot "wave.dat"
