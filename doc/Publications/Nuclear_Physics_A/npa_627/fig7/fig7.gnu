set terminal pslatex
set output "fig7.tex"
set xrange [0.1:1.0]
set yrange [-12.:10.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal U}(n,\chi_p=0)$ (MeV)"
plot "bonn.dat" title "${\cal U}$ with CD-Bonn without $^1S_0$" w line 1, "nim1.dat" title "${\cal U}$ with Nijm I without $^1S_0$" w line 2, "nim2.dat"  title "${\cal U}$ with Nijm II without $^1S_0$" w line 3, "reid.dat"  title "${\cal U}$ with Reid93 without $^1S_0$" w line 4



