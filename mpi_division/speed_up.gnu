set xlabel "nombre de coeur"
set ylabel "Speed_up"
set key box
set key left top
set grid
f(x) = x
g(x) = a*x + b
h(x) = c*x + d
fit g(x) 'speed_up.dat' u  1:3 via a, b
fit h(x) 'speed_up_mpi.dat' u 1:3 via c, d
plot [0:14][0:14]'speed_up.dat' u 1:3 title "Data_OPENMP", f(x) title "f(x) = x", g(x) title "Linear regression", 'speed_up_mpi.dat' u 1:3 title "Data_MPI", h(x) title "Linear regression" 
pause -1
