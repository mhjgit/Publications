set terminal pslatex
set output "model5b.tex"
set xrange [0:6]
set yrange [0.:12]
set xlabel "$\langle E \rangle$"
set ylabel "$n$"
plot "e11b"  title "$\delta=0.5$ for N=11" w l, "e12b"  title "$\delta=0.5$ for N=12" w l



