set terminal pslatex
set output "fig4.tex"
set xrange [8:18]
set yrange [0.4:0.6]
set xlabel "Number of particles $N$"
set ylabel "$\Delta F(E)/N$"
plot "deltaf"  title "$T=0.85$" w p

