set terminal pslatex
set output "fig8c.tex"
set xrange [0.1:1.]
set yrange [0.:25.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$\Delta{\cal S}(n) (MeV)$"
plot "bonn_reid_full" title "$\Delta {\cal S}$ for $l\leq 10$" w line 1, "bonn_reid_s" title "$\Delta{\cal S}$ for $l=0$" w line 2, "bonn_reid_p" title "$\Delta{\cal S}$ for $l\leq 1$" w line 3
