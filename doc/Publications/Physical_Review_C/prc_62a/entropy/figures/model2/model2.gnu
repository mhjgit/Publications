set terminal pslatex
set output "model2.tex"
set pointsize 1
set xrange [0:50]
set yrange [0.:16]
set xlabel "Excitation energy $E$"
set ylabel "Entropy $S/k_B$"
plot "e12a"  title "$\delta=0.1$" w p



