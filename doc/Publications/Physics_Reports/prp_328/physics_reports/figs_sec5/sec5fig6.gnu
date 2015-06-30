set terminal pslatex
set output "sec5fig6.tex"
set xrange [0.3:1.8]
set yrange [9.:20.]
set xlabel "$n_c$ (fm$^{-3}$)"
set ylabel "$R$ km"
plot "star_maxwell100.dat"  title "$B^{1/4}=100$ MeV " w l 1 , "star_maxwell150.dat" title "$B^{1/4}=150$ MeV "  w l 2, "star_maxwell200.dat" title "$B^{1/4}=200$ MeV "  w l 3, "star02.dat" title "$\delta=0.2$"  w l 4






