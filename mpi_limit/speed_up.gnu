set xlabel "Nombres de taches"
set ylabel "Speed_up"
set key box
set key left top
set grid
f(x) = x
plot [0:110][0:110] 'speed_up_domain.d' u 1:3 title "Data_MPI", f(x) title "f(x) = x", 'speed_up_mpi_share.d' u 1:3
pause -1
