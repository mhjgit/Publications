set terminal pslatex
set output "fig4b.tex"
set ylabel "$\left\langle k_F \right| V($^3F_2$)\left| k_F\right\rangle$ (MeVfm$^{-3}$)"
plot "bonnf.dat" w line 1, "nijm1f.dat" w line 2, "nijm2f.dat" w line 3,  "v18f.dat"  w line 4



