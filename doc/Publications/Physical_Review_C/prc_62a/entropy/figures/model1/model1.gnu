set terminal pslatex
set output "model1a.tex"
set xrange [0:12]
set yrange [0.:14]
set xlabel "$E$"
set ylabel "$S/k_B$"
plot "e12a"  title "$\delta=0.5$" w dots, "e12b"  title "$\delta=0.5$ with smoothing" w l



