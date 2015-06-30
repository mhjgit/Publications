set terminal pslatex
set output "fig3.tex"
set xrange [0.3:1.8]
set yrange [0.5:3.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot  "massb150.dat" title "Mixed phase $B^{1/4}=150$ MeV "  w l 1, "massb200.dat" title "Mixed phase $B^{1/4}=200$ MeV "  w l 2, "mass_maxwell150.dat" title "Maxwell $B^{1/4}=150$ MeV "  w l 3, "mass_maxwell200.dat" title "Maxwell $B^{1/4}=200$ MeV "  w l 4, "rotm200.dat"  title "Rotational mass for $B^{1/4}=200$ MeV " w l 5 






