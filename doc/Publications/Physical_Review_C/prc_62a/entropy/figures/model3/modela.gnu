set terminal pslatex
set output "model3a.tex"
set xrange [0:5]
set yrange [0.:18]
set xlabel "$k_BT$"
set ylabel "$S/k_B$"
plot "s11a"  title "$\delta=0.1$ for N=11" w l, "s12a"  title "$\delta=0.1$ for N=12" w l



