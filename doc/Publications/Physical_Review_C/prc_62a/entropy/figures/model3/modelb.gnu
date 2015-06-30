set terminal pslatex
set output "model3b.tex"
set xrange [0:1.5]
set yrange [0.:18]
set xlabel "$k_BT$"
set ylabel "$S/k_B$"
plot "s11b"  title "$\delta=0.5$ for N=11" w l, "s12b"  title "$\delta=0.5$ for N=12" w l



