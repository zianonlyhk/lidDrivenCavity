#! /opt/homebrew/bin/gnuplot

cd '/Users/zianhuang/Room214N/dev/mphil/lidDrivenCavity/data'

set terminal gif animate delay 5
set output 'gif/uPlot.gif'
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_uResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_uResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}

reset session

cd '/Users/zianhuang/Room214N/dev/mphil/lidDrivenCavity/data'

set terminal gif animate delay 5
set output 'gif/vPlot.gif'
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_vResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_vResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}

reset session

cd '/Users/zianhuang/Room214N/dev/mphil/lidDrivenCavity/data'

set terminal gif animate delay 5
set output 'gif/pPlot.gif'
set xyplane at 0
set xlabel "x"
set ylabel "y"
stats 'testSolver_pResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot 'testSolver_pResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}