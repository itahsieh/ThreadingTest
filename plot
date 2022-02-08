#!/usr/bin/gnuplot

DatFile='ParallelEfficiency.dat'

set terminal png
set output 'efficiency.png'

set xlabel 'Number of Threads'
set ylabel 'Efficiency'
set xrange[1:10]
set yrange[0:1]
set linetype 1 lc rgb "red" lw 2 ps 3 pt 7
set linetype 2 lc rgb "blue" lw 2 ps 3
plot DatFile index 1 using 1:3 lt 1 smooth acsplines title 'OpenMP',\
     DatFile index 1 using 1:3 lt 1 notitle,\
     DatFile index 2 using 1:3 lt 2 smooth acsplines title 'PThread',\
     DatFile index 2 using 1:3 lt 2 notitle

set output 'speedup.png'
set ylabel 'Speed-up'
set yrange[1:4]
plot DatFile index 1 using 1:5 lt 1 smooth acsplines title 'OpenMP',\
     DatFile index 1 using 1:5 lt 1 notitle,\
     DatFile index 2 using 1:5 lt 2 smooth acsplines title 'PThread',\
     DatFile index 2 using 1:5 lt 2 notitle
