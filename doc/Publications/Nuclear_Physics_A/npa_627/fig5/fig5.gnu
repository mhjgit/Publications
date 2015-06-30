set terminal pslatex
set output "fig5.tex"
set xrange [0.1:1.0]
set yrange [-5.:10.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "${\cal U}(n,\chi_p=0.5)$ (MeV)"
plot "bonn.dat" title "${\cal U}$ with CD-Bonn for $l\geq 1$" w line 1, "nim1.dat" title "${\cal U}$ with Nijm I for $l\geq 1$" w line 2, "nim2.dat"  title "${\cal U}$ with Nijm II for $l\geq 1$" w line 3, "reid.dat"  title "${\cal U}$ with Reid93 for $l\geq 1$ " w line 4, "v18.dat"  title "${\cal U}$ with Argonne$V_{18}$ for $l\geq 1$ " w line 5



