set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
set size ratio -1
set cbrange[-0.1:0.1]
set pm3d map
set xlabel "x1"
set ylabel "x2"
set cblabel "sound pressure in time domain"
unset key
set term gif animate delay 0.01
set output "movie.gif"
do for [i = 0:313 ] {
   t=i*0.01
   w=3.0
   splot "domain.res" u ($1):($2):($3)*cos(w*t)+($4)*sin(w*t)
}