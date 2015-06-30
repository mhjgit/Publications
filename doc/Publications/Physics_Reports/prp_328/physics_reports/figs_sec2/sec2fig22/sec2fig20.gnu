set terminal pslatex
set output "sec2fig20.tex"
set xrange [0.1:2.0]
set yrange [0.:3.]
set xlabel "$n$ (fm$^{-3}$)"
set ylabel "$(v_s/c)^2$"
plot "inapr" title "Argonne $V_{18}+\delta v +$UIX$\ast$ " w l 1, "in00" title "$\delta=0.0$" w l 2, "in01" title "$\delta=0.1$" w l 3, "in02" title "$\delta=0.2$" w l 4, "in03" title "$\delta=0.3$" w l 5,  "inbk" title "Baym and Kalogera" w l 6



