set terminal pslatex
set output "model5c.tex"
set xrange [0:0.2]
set yrange [0.:14]
set xlabel "$\langle E \rangle$"
set ylabel "$n$"
plot "e11c"  title "$\delta=10$ for N=11" w l, "e12c"  title "$\delta=10$ for N=12" w l



