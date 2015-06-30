set terminal pslatex
set output "sec5fig3.tex"
set xrange [0.3:1.8]
set yrange [0.5:3.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "massb100.dat"  title "$B^{1/4}=100$ MeV " w l 1 , "massb150.dat" title "$B^{1/4}=150$ MeV "  w l 2, "massb200.dat" title "$B^{1/4}=200$ MeV "  w l 3, "mass02.dat" title "$\delta=0.2$"  w l 4



