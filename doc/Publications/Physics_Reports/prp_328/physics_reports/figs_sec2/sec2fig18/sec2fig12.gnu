set terminal pslatex
set output "sec2fig15.tex"
set xrange [0.1:0.7]
set yrange [0.:250]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "PNM ${\cal E}$ (MeV)"
plot "v18dv_pnm.dat" title "Argonne $V_{18}+\delta v$ " w l 1, "v18dvuix_pnm.dat" title "Argonne $V_{18}+\delta v+UIX$" w l 2, "bhf_pnm.dat" title "BHF" w l 3, "dbhf_pnm.dat" title "DBHF" w l 4
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set yrange [-25.:150.]
set ylabel "SNM ${\cal E}$ (MeV)"
plot "v18dv_snm.dat" title "Argonne $V_{18}+\delta v$ " w l 1, "v18dvuix_snm.dat" title "Argonne $V_{18}+\delta v+UIX$" w l 2, "bhf_snm.dat" title "BHF" w l 3, "dbhf_snm.dat" title "DBHF" w l 4



