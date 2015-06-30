set terminal pslatex
set output "fig1a.tex"
set xrange [0.:1.2]
set yrange [800:1300]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\mu_i$ (MeV)"
plot "bhf_n.dat" title "$\mu_n$" w l 1, "bhf_p.dat" title "$\mu_p$" w l 2, "bhf_sm.dat" title "$\mu_{\Sigma^-}$" w l 3, "bhf_l.dat" title "$\mu_{\Lambda}$" w l 4, "bhf_s0f.dat" title "$\mu_{\Sigma^0}$" w l 5, "bhf_spf.dat" title "$\mu_{\Sigma^{+}}$" w l 6



