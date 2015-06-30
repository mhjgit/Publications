set terminal pslatex
set output "fig3a.tex"
set logscale y
set xrange [0.:1.2]
set ylabel "Fraction"
plot "bhf_n.dat" title "$n$" w l 1, "bhf_p.dat" title "$p$" w l 2, "bhf_e.dat" title "$e^-$" w l 3, "bhf_mu.dat" title "$\mu^-$" w l 4, "bhf_sm.dat" title "$\Sigma^-$" w l 5
