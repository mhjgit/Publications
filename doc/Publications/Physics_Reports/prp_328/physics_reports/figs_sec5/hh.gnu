set terminal pslatex
set output "hh.tex"
set xrange [9.:17.]
set yrange [0.5:2.5]
set xlabel "$R$ km"
set ylabel "$M/M_{\odot}$"
plot "rm02.dat" title "$\delta=0.2$" w l 1, "rm150.dat" title "$B^{1/4}=150$ MeV "  w l 2, "rm200.dat" title "$B^{1/4}=200$ MeV "  w l 3






