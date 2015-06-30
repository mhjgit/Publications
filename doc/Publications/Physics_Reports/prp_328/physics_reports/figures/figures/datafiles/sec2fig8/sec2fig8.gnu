set terminal pslatex
set output "sec2fig8.tex"
set xrange [0.1:1.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "bonn_neut.dat" title "PNM" w l 1, "bonn_beta.dat" title "$\beta$-stable" w l 2
set xrange [0.1:1.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$P$ (MeVfm$^{-3}$)"
plot "pneut.dat" title "PNM" w l 1, "pbeta.dat" title "$\beta$-stable" w l 2

