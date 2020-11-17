# gnuplot -e "datafile='contr_tot_gdr12.dat'; outfile='fig_pdf12.pdf'; nneigh=6" plot_contr_pdf.gnuplot

set terminal pdf enhanced font ",12"
set output outfile

set xlabel "Neighbor index"
set ylabel "{/Symbol s}"

unset key

plot datafile u 1:2 pt 7 ps 0.5

#plot datafile u 1:2 w l lw 2 title "total", \
#     for [t = 3:2+nneigh] datafile u 1:t w l title "pair ".(t-2)