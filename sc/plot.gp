set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
set size ratio -1
set cbrange[-0.1:0.1]
set pm3d map
set term png
set output "result.png"
set xlabel "x1"
set ylabel "x2"
set cblabel "real part of complex sound pressure"
sp "domain.res"