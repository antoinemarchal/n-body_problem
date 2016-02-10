set xlabel "iteration"
set ylabel "Energie"
set key box
set key left top
plot 'output.dat' u 1:2 w l title "V", 'output.dat' u 1:3 w l title "U", 'output.dat' u 1:4 w l title "E_{tot}" 
pause -2
