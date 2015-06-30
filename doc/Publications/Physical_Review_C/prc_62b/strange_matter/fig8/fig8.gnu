set terminal pslatex
set output "fig8.tex"
set xrange [0.3:2.]
set yrange [10.:120]
set xlabel "$M/M_{\odot}$"
set ylabel "$I$"
plot "i02.dat" title "APR98 pn-matter"  w l 1, "ihyp.dat" title "APR98 with hyperons"  w l 2
