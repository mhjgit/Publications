set terminal pslatex
set output "fig4c.tex"
set xlabel "$k_F$ (fm$^{-1}$)"
set ylabel "$\left\langle k_F \right| V({^}3P_2--{^}3F_2)\left| k_F\right\rangle$ (MeVfm$^{-3}$)"
plot "bonnpf.dat" w line 1, "nijm1pf.dat" w line 2, "nijm2pf.dat" w line 3,  "v18pf.dat" w line 4



