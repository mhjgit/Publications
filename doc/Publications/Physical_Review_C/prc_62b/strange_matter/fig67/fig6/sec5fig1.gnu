set terminal pslatex
set output "fig5.tex"
set xrange [0.2:1.8]
set yrange [0.5:2.5]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "mass02_nr"  title "pn-interaction" w l 1 , "mass_allYY_nr" title "YY-interaction$\Delta M_r$"  w l 2,"mass_noYY_nr" title "YN-interaction"  w l 3



