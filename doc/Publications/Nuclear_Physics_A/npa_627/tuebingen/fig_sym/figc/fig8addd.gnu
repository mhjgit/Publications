set terminal pslatex
set output "fig8addd.tex"
set xrange [0.1:1.]
set yrange [0.:70.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal S}(n) (MeV)$"
plot "bonn_p" title "${\cal S}$ with CD-Bonn" w line 1, "nim1_p" title "${\cal S}$ with Nijm I" w line 2, "nim2_p"  title "${\cal S}$ with Nijm II" w line 3,  "reid_p"  title "${\cal S}$ with Reid93" w line 4



