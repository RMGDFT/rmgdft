#!/usr/bin/gnuplot

set encoding iso_8859_1

#set terminal postscript {landscape | portrait | eps | default}
#                        {enhanced | noenhanced}
#                        {color | monochrome} {solid | dashed}
#                        {<duplexing>}
#                        {"<fontname>"} {<fontsize>}
#set terminal postscript eps enhanced color solid 'Times-Roman' 18
#set output 'test.edos.eps'
set term png
set output 'edos.png'



#set border 31 lt 7 lw 2
#set zeroaxis lt 7 lw 1
#set grid     lt 7 lw 1

#set size 1.0, 1.3

#set lmargin  9; set rmargin  9
#set tmargin  2; set bmargin  4



#set title ''    font 'Times-Roman, 15'
#set label 2 '(c)' at graph l2x, l2y right font 'Times-Roman, 18'
# Angstrom -> \305 

EF = 0

Emin = EF-2
Emax = EF+2
dE   = 0.4
Etics = 5

xmin = 0
xmax = 656

ymin = 0
ymax = 80

set samples 25


set pm3d map

#          * there are 37 available rgb color mapping formulae:
#             0: 0               1: 0.5             2: 1
#             3: x               4: x^2             5: x^3
#             6: x^4             7: sqrt(x)         8: sqrt(sqrt(x))
#             9: sin(90x)       10: cos(90x)       11: |x-0.5|
#            12: (2x-1)^2       13: sin(180x)      14: |cos(180x)|
#            15: sin(360x)      16: cos(360x)      17: |sin(360x)|
#            18: |cos(360x)|    19: |sin(720x)|    20: |cos(720x)|
#            21: 3x             22: 3x-1           23: 3x-2
#            24: |3x-1|         25: |3x-2|         26: (3x-1)/2
#            27: (3x-2)/2       28: |(3x-1)/2|     29: |(3x-2)/2|
#            30: x/0.32-0.78125 31: 2*x-0.84       32: 4x;1;-2x+1.84;x/0.08-11.5
#            33: |2*x - 0.5|    34: 2*x            35: 2*x - 0.5
#            36: 2*x - 1
#          * negative numbers mean inverted=negative colour component
#          * thus the ranges in `set pm3d rgbformulae' are -36..36


#    7,5,15   ... traditional pm3d (black-blue-red-yellow)
#    3,11,6   ... green-red-violet
#    23,28,3  ... ocean (green-blue-white); try also all other permutations
#    21,22,23 ... hot (black-red-yellow-white)
#    30,31,32 ... color printable on gray (black-blue-violet-yellow-white)
#    33,13,10 ... rainbow (blue-green-yellow-red)
#    34,35,36 ... AFM hot (black-red-yellow-white)

#set palette rgbformulae 21,22,23
set palette rgbformulae 33,13,10
#set palette rgbformulae 7,5,15
#set palette rgbformulae 34,35,36
#set palette functions sin(0.5*pi*gray), gray**4, 0*gray

set xtics 

set ytics 1

#set key top left
set nokey

#set ytics dE
#set ytics nomirror 
#set mytics Etics
#set noxtics 
#set mxtics 0
#set noy2tics

#set format cb '10^{%.0f}'
set format cb '%.1f'
set cbrange [0:100.0]

#set cbrange [0:1]

#set yrange [-4:0]
#set format y '%.0f'
#set ylabel 'Energy [eV]' font 'Times-Roman, 20'

#set xrange [-15:15]
#set format x '%.0f'
#set xlabel 'x [\305]' font 'Times-Roman, 20'


#splot 'test.edos' u ($1):($2-EF):(log($3)/log(10))
splot './dos.dat' u ($1):($2-EF):($3)
#splot [] [] 'dos.dat' u (log($1))



     


     #     'mbp.scf-01.eig' u (dx+$3/ne*2*dx+$1):4:(0.1) not w yerrorbars lt 1 lw 1 pt 0,\

#guide for line and point styles:
#  0  ..............  .                    broken line
#  1  --------------  +                    red
#  2  -- -- -- -- --  x                    green
#  3  -  -  -  -  -   *                    blue
#  4  ..............  empty square         magenta
#  5  __.__.__.__.__  full  square         cyan
#  6  _ . _ . _ . _   empty circle         yellow
#  7  - -  - -  - -   full  circle         black
#  8  - - -  - - -    empty up triangle    brown
#  9  - - - -  - - -  full  up triangle    grey
# 10 (1)              empty down triangle
# 11 (2)              full  down triangle
# 12 (3)              empty diamond
# 13 (4)              full  diamond
# 14 (5)              empty pentagon
# 15 (6)              full  pentagon
# 16-31               watches



# ------------
#    done
# ------------



