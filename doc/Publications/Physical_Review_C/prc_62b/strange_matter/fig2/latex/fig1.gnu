set terminal pslatex
set output "fig2a.tex"
set xrange [0.2:1.2]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$U_i^N$ (MeV)"
plot "chemf_nn.dat" title "$U^N_n$" w l 1, "chemf_pn.dat" title "$U^N_p$" w l 2, "chemf_smn.dat" title "$U^N_{\Sigma^-}$" w l 3, "chemf_ln.dat" title "$U^N_{\Lambda}$" w l 4



