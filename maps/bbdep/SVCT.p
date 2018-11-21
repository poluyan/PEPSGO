reset
set terminal pngcairo size 1024, 768 enhanced font 'Verdana,12'
set output 'SVCT.png'

set style line 11 lc rgb '#000000' lt 2 lw 1
set border 1+2+4+8 ls 11

# define grid
set style line 12 lc rgb '#e1e1e1' lt 1 lw 1
set grid back ls 12

# color definitions
set style line 1 lt rgb 'red' lw 1.0 ps 1
set style line 1 lt rgb 'green' lw 1.0 ps 1
set style line 1 lt rgb 'blue' lw 1.0 ps 1
set style line 1 lt rgb 'black' lw 1.0 ps 1

set xrange [0:360]

plot \
'SVCT.dat' using 1:2 with lines ls 1 title 'ser',\
'SVCT.dat' using 1:3 with lines ls 2 title 'val',\
'SVCT.dat' using 1:4 with lines ls 3 title 'cys',\
'SVCT.dat' using 1:5 with lines ls 4 title 'thr'