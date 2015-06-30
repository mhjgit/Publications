set terminal pslatex
set output "sec2fig12.tex"
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "argneut.dat" title "Argonne $V_{18}$ LOB" w l 1, "v18uix_pnm.dat" title "Argonne $V_{18}$ VCS+UIX" w l 2, "baldo_pnm.dat" title "3 hole-line" w l 3
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "argsnm.dat" title "Argonne $V_{18}$ LOB" w l 1, "v18uix_snm.dat" title "Argonne $V_{18}$ VCS+UIX" w l 2, "baldo_snm.dat" title "3 hole-line" w l 3



