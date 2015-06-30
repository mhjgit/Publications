set terminal pslatex
set output "sec5fig3.tex"
set xrange [0.5:3.]
set yrange [0.:120]
set xlabel "$M/M_{\odot}$"
set ylabel "$I$"
plot "i013.dat"  title "$\delta=0.13$" w l 1 , "i02.dat" title "$\delta=0.2$"  w l 2, "i03.dat" title "$\delta=0.3$"  w l 3, "i04.dat" title "$\delta=0.4$"  w l 4



