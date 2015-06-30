set terminal pslatex
set output "fig8.tex"
set xrange [0.1:1.]
set yrange [0.:100]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal S}(n)$ (MeV)"
plot "bonn_sym.dat" title "${\cal S}$ with CD-Bonn" w line 1, "nim1_sym.dat" title "${\cal S}$ with Nijm I" w line 2, "nim2_sym.dat"  title "${\cal S}$ with Nijm II" w line 3,  "reid_sym.dat"  title "${\cal S}$ with Reid93" w line 4,  "argonne_sym.dat"  title "${\cal S}$ with $V_{18}$" w line 5



