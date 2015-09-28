set terminal pslatex
set output "sec2fig9.tex"
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$v_s/c$"
plot "vcneut.dat" title "PNM" w l 1, "vcbeta.dat" title "$\beta$-stable" w l 2
set xrange [0.1:0.7]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\Gamma$"
plot "gammaneut.dat" title "PNM" w l 1, "gammabeta.dat" title "$\beta$-stable" w l 2