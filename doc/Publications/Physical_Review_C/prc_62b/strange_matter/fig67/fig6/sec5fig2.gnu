set terminal pslatex
set output "fig7.tex"
set xrange [9.:14.]
set yrange [0.5:2.1]
set xlabel "$R$ km"
set ylabel "$M/M_{\odot}$"
plot "mr02_nr"  title "pn-interaction" w l 1 , "mr_all.dat" title "YY-interaction"  w l 2, "mr_no.dat" title "YN-interaction"  w l 3



