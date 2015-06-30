set terminal pslatex
set output "model5a.tex"
set xrange [0:30]
set yrange [0.:12]
set xlabel "$\langle E \rangle$"
set ylabel "$n$"
plot "e11a"  title "$\delta=0.1$ for N=11" w l, "e12a"  title "$\delta=0.1$ for N=12" w l



