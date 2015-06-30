set terminal pslatex
set output "sec5fig9.tex"
set xrange [0.3:1.8]
set yrange [0.:150]
set xlabel "$n_c$ fm$^{-3}$"
set ylabel "$I$"
plot "di02.dat" title "$\delta=0.2$"  w l 1, "di200.dat" title "Mixed phase $B^{1/4}=200$ MeVfm$^{-3}$ "  w l 2, "dimaxwell200.dat" title "Maxwell $B^{1/4}=200$ MeVfm$^{-3}$ "  w l 3 
