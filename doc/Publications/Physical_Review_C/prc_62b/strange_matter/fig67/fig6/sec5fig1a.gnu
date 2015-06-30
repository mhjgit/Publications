set terminal pslatex
set output "fig6.tex"
set xrange [0.2:1.8]
set yrange [0.5:3.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "mass02_r" title "pn-interaction"  w l 1, "mass_allYY_r" title "YY-interaction"  w l 2, "mass_noYY_r" title "YN-interaction"  w l 3



