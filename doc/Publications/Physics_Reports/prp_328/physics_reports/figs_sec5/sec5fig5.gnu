set terminal pslatex
set output "sec5fig5.tex"
set xrange [0.3:1.8]
set yrange [0.5:3.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "mass_maxwell100.dat"  title "$B^{1/4}=100$ MeVfm$^{-3}$ " w l 1 , "mass_maxwell150.dat" title "$B^{1/4}=150$ MeVfm$^{-3}$ "  w l 2, "mass_maxwell200.dat" title "$B^{1/4}=200$ MeVfm$^{-3}$ "  w l 3



