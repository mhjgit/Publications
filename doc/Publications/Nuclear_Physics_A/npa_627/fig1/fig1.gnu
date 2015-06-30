set terminal pslatex
set output "fig1.tex"
set xrange [0.1:0.72]
set yrange [10.:80]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal S}(n) (MeV)$"
plot "paris.dat" title "${\cal S}$ with Paris" w line 1, "reid.dat" title "${\cal S}$ with Reid68" w line 2, "v14.dat"  title "${\cal S}$ with $V_{14}$" w line 3



