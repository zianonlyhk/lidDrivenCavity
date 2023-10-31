#! /opt/homebrew/bin/gnuplot

x0 = 0
x1 = 1
y0 = 0
y1 = 1

cd './data'
set size square
set terminal gif animate delay 15
set output 'gif/uPlot.gif'
set title "x-velocity"
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_uResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set xrange [x0:x1]
set yrange [y0:y1]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_uResults.dat' index (i-1) using 2:3:4 with image palette notitle
}

reset 

set size square
set terminal gif animate delay 15
set output 'gif/vPlot.gif'
set title "y-velocity"
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_vResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set xrange [x0:x1]
set yrange [y0:y1]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_vResults.dat' index (i-1) using 2:3:4 with image palette notitle
}

reset 

set size square
set terminal gif animate delay 15
set output 'gif/pPlot.gif'
set title "pressure"
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_pResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set xrange [x0:x1]
set yrange [y0:y1]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_pResults.dat' index (i-1) using 2:3:4 with image palette notitle
}

# reset 

# set size square
# set terminal gif animate delay 15
# set output 'gif/arrowPlot.gif'
# set xyplane at 0
# set xlabel "x"
# set ylabel "y"
# set xrange [x0:x1]
# set yrange [y0:y1]
# stats 'testSolver_vecResults.dat' using 4 nooutput
# do for [i=1:int(STATS_blocks)-1] {
#     plot "testSolver_vecResults.dat" index (i-1) using 2:3:4:5 with vectors filled head lw 0.75 notitle
# }

# reset session

# set size square
# set terminal gif animate delay 15
# set output 'gif/streamFuncPlot.gif'
# set xyplane at 0
# set xlabel "x"
# set ylabel "y"
# stats 'testSolver_streamFuncResults.dat' using 4 nooutput
# set cbrange [STATS_min:STATS_max]
# set xrange [x0:x1]
# set yrange [y0:y1]
# do for [i=1:int(STATS_blocks)-1] {
#     plot 'testSolver_streamFuncResults.dat' index (i-1) using 2:3:4 with image palette notitle
# }

reset 

set size square
set terminal gif animate delay 15
set output 'gif/arrAndSfPlot.gif'
set title "velocity direction and magnitude"
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_streamFuncResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set xrange [x0:x1]
set yrange [y0:y1]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_streamFuncResults.dat' index (i-1) using 2:3:4 with image palette notitle,\
    "testSolver_vecResults.dat" index (i-1) using 2:3:4:5 with vectors filled head lw 0.8 lc rgb 0x056608 notitle
}

reset 

set size square
set terminal gif animate delay 15
set output 'gif/tempPlot.gif'
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_tempResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set xrange [x0:x1]
set yrange [y0:y1]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_tempResults.dat' index (i-1) using 2:3:4 with image palette notitle
}
