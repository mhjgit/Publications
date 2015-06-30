set terminal pslatex
set output "fig4.tex"
set xrange [0.1:1.0]
set yrange [-40.:-10.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal U}(n,\chi_p=0.5)$ (MeV)"
plot "bonn.dat" title "${\cal U}$ with CD-Bonn for $^1S_0$" w line 1, "nim1.dat" title "${\cal U}$ with Nijm I for$^1S_0$" w line 2, "nim2.dat"  title "${\cal U}$ with Nijm II for $^1S_0$" w line 3, "reid.dat"  title "${\cal U}$ with Reid93 for $^1S_0$" w line 4, "v18.dat"  title "${\cal U}$ with Argonne $V_{18}$ for $^1S_0$" w line 5



