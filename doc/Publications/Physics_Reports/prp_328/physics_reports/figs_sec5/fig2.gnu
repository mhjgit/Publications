set terminal pslatex
set output "fig2.tex"
set xrange [0.3:1.8]
set yrange [0.5:3.2]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$M/M_{\odot}$"
plot "mass013.dat"  title "$\delta=0.13$" w l 1 , "mass02.dat" title "$\delta=0.2$"  w l 2, "mass03.dat" title "$\delta=0.3$"  w l 3, "mass04.dat" title "$\delta=0.4$"  w l 4, "rotm02.dat" title "Rotational mass for $\delta=0.2$"  w l 5





