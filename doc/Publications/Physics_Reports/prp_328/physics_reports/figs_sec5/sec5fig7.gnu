set terminal pslatex
set output "sec5fig7.tex"
set xrange [0.4:1.8]
set yrange [0.5:3.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "rotm200.dat"  title "Rotational mass for $B^{1/4}=200$ MeV " w l 1 , "massb200.dat" title "$B^{1/4}=200$ MeV "  w l 2, "mass02.dat" title "$\delta=0.2$"  w l 3, "rotm02.dat" title "Rotational mass for $\delta=0.2$"  w l 4
