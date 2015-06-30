set terminal pslatex
set output "sec2fig18.tex"
set xrange [0.1:0.8]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "PNM ${\cal E}$ (MeV)"
plot "v18dvuix_pnm.dat" title "Argonne $V_{18}+\delta v+UIX$" w l 1, "param_pnm.dat" title "Parametrization" w l 2
set xrange [0.1:0.8]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "SNM ${\cal E}$ (MeV)"
plot "v18dvuix_snm.dat" title "Argonne $V_{18}+\delta v+UIX$" w l 1, "param_snm.dat" title "Parametrization" w l 2




