set terminal pslatex
set output "model3c.tex"
set xrange [0:1]
set yrange [0.:18]
set xlabel "$k_BT$"
set ylabel "$S/k_B$"
plot "s11c"  title "$\delta=10$ for N=11" w l, "s12c"  title "$\delta=10$ for N=12" w l



