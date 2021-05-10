#CONTOUR B(R,s) PLOT

reset
unset key

set xlabel "R"
set ylabel "s"
set xrange [20:224]
set yrange [3:59]
#set zrange [-200:600]
#set title "B(R,s)"

#set pm3d map
set view map
#set dgrid3d
set contour
unset surface
#set palette rgbformulae 33,13,10
#set cbrange[-200:200]
set cntrparam levels incremental -2,0.05,2 #for 4
set cbrange[-2:2]
#set cntrparam levels incremental -40,1,0 #for 4

splot "Brs.dat" u 1:2:3 palette notitle w l
