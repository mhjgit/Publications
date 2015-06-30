set terminal pslatex
set output "model4.tex"
set xrange [0:25]
set yrange [0.:5]
set xlabel "$E$"
set ylabel "$\Delta S/k_B$"
plot "s005"  title "$\delta=0.05$" w l, "s1"  title "$\delta=1$" w l, "s5"  title "$\delta=5$" w l



