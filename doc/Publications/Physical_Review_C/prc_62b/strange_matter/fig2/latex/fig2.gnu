set terminal pslatex
set output "fig2b.tex"
set xrange [0.2:1.2]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$U_i^Y$ (MeV)"
plot "chemf_ns.dat" title "$U^Y_n$" w l 1, "chemf_ps.dat" title "$U^Y_p$" w l 2, "chemf_sms.dat" title "$U^Y_{\Sigma^-}$" w l 3, "chemf_ls.dat" title "$U^Y_{\Lambda}$" w l 4



