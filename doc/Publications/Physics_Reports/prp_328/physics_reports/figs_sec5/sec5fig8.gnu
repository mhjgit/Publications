set terminal pslatex
set output "sec5fig8.tex"
set xrange [1.2:2.]
set yrange [40.:102]
set xlabel "$M/M_{\odot}$"
set ylabel "$I$"
plot "i02.dat" title "$\delta=0.2$"  w l 1, "ib200.dat" title "$B^{1/4}=200$ MeV "  w l 2, "ib150.dat" title "$B^{1/4}=150$ MeV "  w l 3, "ib100.dat" title "$B^{1/4}=100$ MeV "  w l 4
