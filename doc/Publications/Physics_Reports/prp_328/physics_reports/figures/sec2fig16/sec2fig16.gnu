set terminal pslatex
set output "sec2fig16.tex"
set xrange [0.1:1.0]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "apr_enuc.dat" title "Argonne $V_{18}+\delta v +UIX\ast$ " w l 1, "dbhf_enuc.dat" title "DBHF" w l 2
set xrange [0.1:1.0]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\chi_p$"
plot "apr_prot.dat" title "Argonne $V_{18}+\delta v+UIX\ast$ " w l 1, "dbhf_prot.dat" title "DBHF" w l 2,"bhf_prot.dat" title "BHF" w l 3, "el.dat" title "$e^-$ direct URCA" w l 4, "mu.dat" title "$\mu^-$ direct URCA" w l 5



