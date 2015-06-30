set terminal pslatex
set output "fig4a.tex"
set ylabel "$\left\langle k_F \right| V(^3P_2)\left| k_F\right\rangle$ (MeVfm$^{-3}$)"
plot "bonndiag.dat" title "CD-Bonn" w line 1, "nijm1diag.dat" title "Nijm I" w line 2, "nijm2diag.dat"  title "Nijm II" w line 3,  "v18diag.dat"  title "Argonne $V_{18}$" w line 4



