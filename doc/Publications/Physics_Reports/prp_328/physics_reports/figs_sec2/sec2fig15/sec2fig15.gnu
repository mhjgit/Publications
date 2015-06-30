set terminal pslatex
set output "sec2fig15.tex"
set xrange [1.0:2.0]
set xlabel "$k_F$ (fm$^{-1}$)"
set ylabel "${\cal E}$ (MeV)"
plot "bub.dat" title "BUB" w l 1, "ring.dat" title "RING" w l 2, "high.dat" title "HIGH" w l 3, "total.dat" title "TOTAL" w l 4

