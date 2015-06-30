set terminal pslatex
set output "sec2fig19.tex"
set xrange [0.1:1.0]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "apr_enuc.dat" title "Argonne $V_{18}+\delta v +UIX\ast$ " w l 1, "enuc00.dat" title "$\delta=0.0$" w l 2, "enuc01.dat" title "$\delta=0.1$" w l 3, "enuc02.dat" title "$\delta=0.2$" w l 4, "enuc03.dat" title "$\delta=0.3$" w l 5
set xrange [0.1:1.0]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\chi_p$"
plot "apr_prot.dat" title "Argonne $V_{18}+\delta v +UIX\ast$ " w l 1, "prot00.dat" title "$\delta=0.0$" w l 2, "prot01.dat" title "$\delta=0.1$" w l 3, "prot02.dat" title "$\delta=0.2$" w l 4, "prot03.dat" title "$\delta=0.3$" w l 5



