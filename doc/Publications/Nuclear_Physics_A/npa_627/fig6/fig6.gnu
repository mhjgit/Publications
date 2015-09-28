set terminal pslatex
set output "fig6.tex"
set xrange [0.1:1.0]
set yrange [-33.:-1.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal U}(n,\chi_p=0)$ (MeV)"
plot "bonn.dat" title "${\cal U}$ with CD-Bonn" w line 1, "nim1.dat" title "${\cal U}$ with Nijm I" w line 2, "nim2.dat"  title "${\cal U}$ with Nijm II" w line 3, "reid.dat"  title "${\cal U}$ with Reid93" w line 4, "v18.dat"  title "${\cal U}$ with Argonne $V_{18}$" w line 5


