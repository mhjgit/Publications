set terminal pslatex
set output "fig1b.tex"
set xrange [0.:1.2]
set yrange [800:1300]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\mu_i$ (MeV)"
plot "apr_n.dat" title "$\mu_n$" w l 1, "apr_p.dat" title "$\mu_p$" w l 2, "apr_sm.dat" title "$\mu_{\Sigma^-}$" w l 3, "apr_l.dat" title "$\mu_{\Lambda}$" w l 4, "apr_s0f.dat" title "$\mu_{\Sigma^0}$" w l 5, "apr_spf.dat" title "$\mu_{\Sigma^{+}}$" w l 6



