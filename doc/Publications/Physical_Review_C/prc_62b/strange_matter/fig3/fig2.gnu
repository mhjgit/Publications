set terminal pslatex
set output "fig3b.tex"
set logscale y
set xrange [0.:1.2]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "Fraction"
plot "apr_n.dat" title "$n$" w l 1, "apr_p.dat" title "$p$" w l 2, "apr_e.dat" title "$e^-$" w l 3, "apr_mu.dat" title "$\mu^-$" w l 4, "apr_sm.dat" title "$\Sigma^-$" w l 5, "apr_l.dat" title "$\Lambda$" w l 6



