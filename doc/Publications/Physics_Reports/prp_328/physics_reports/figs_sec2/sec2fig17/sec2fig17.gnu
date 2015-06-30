set terminal pslatex
set output "sec2fig17.tex"
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "PNM ${\cal E}$ (MeV)"
plot "v18.dat" title "Argonne $V_{18}$ VCS" w l 1, "v18dvuix_pnm.dat" title "Argonne $V_{18}+$UIX" w l 2, "baldo_pnm" title "3-hole line" w l 3
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "SNM ${\cal E}$ (MeV)"
plot "v18snm.dat" title "Argonne $V_{18}$ VCS" w l 1, "v18dvuix_snm.dat" title "Argonne $V_{18}+$UIX" w l 2, "baldo_snm" title "3-hole line" w l 3



