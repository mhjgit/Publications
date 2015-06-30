set terminal pslatex
set output "sec3fig1.tex"
set xrange [0.1:1.8]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$P$ MeVfm$^{-3}$"
plot "p200" title "Mixed phase" w l 1, "p200_max" title "Maxwell construction" w l 2

