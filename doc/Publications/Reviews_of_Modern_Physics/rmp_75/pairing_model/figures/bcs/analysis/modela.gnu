set terminal pslatex
set output "fig1.tex"
set xrange [0:25]
set xlabel "Excitation energy $E$"
set ylabel "Free energy $F(E)$"
plot "in50"  title "$T=0.5$" w p
set output "fig2.tex"
set xrange [0:25]
set xlabel "Excitation energy $E$"
set ylabel "Free energy $F(E)$"
plot "in100"  title "$T=1.0$" w p
set output "fig3.tex"
set xrange [0:25]
set xlabel "Excitation energy $E$"
set ylabel "Free energy $F(E)$"
plot "in85"  title "$T=0.85$" w p



