set terminal pslatex
set output "model1b.tex"
set xrange [0:11]
set yrange [0:14]
set xlabel "$E$"
set ylabel "$S/k_B$"
plot "e11a"  title "$\delta=0.5$" w dots, "e11b"  title "$\delta=0.5$ with smoothing" w l


