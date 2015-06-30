set terminal pslatex
set output "fig4.tex"
set xrange [0.:1.2]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal E}$ (MeV)"
plot "bhfpn.dat" title "BHF pn-matter" w l 1, "bhfhyp.dat" title "BHF with hyperons" w l 2, "aprpn.dat" title "APR98 pn-matter" w l 3, "aprhyp.dat" title "APR98 with hyperons" w l 4



