from pymol.cgo import *
from pymol import cmd


def drawbox(x0, y0, z0, x1, y1, z1, linewidth=2.0, r=1.0, g=1.0, b =1.0):
    box = [
        LINEWIDTH, float(linewidth),

        BEGIN, LINES,
        COLOR, float(r), float(g), float(b),

        VERTEX, x0, y0, z0,
        VERTEX, x0, y0, z1,

        VERTEX, x0, y0, z0,
        VERTEX, x0, y1, z0,

        VERTEX, x0, y0, z0,
        VERTEX, x1, y0, z0,

        VERTEX, x1, y1, z0,
        VERTEX, x1, y0, z0,

        VERTEX, x1, y1, z0,
        VERTEX, x0, y1, z0,

        VERTEX, x1, y1, z0,
        VERTEX, x1, y1, z1,

        VERTEX, x1, y0, z1,
        VERTEX, x1, y1, z1,

        VERTEX, x1, y0, z1,
        VERTEX, x1, y0, z0,
         
        VERTEX, x1, y0, z1,
        VERTEX, x0, y0, z1,
         
        VERTEX, x0, y1, z1,
        VERTEX, x0, y0, z1,

        VERTEX, x0, y1, z1,
        VERTEX, x1, y1, z1,

        VERTEX, x0, y1, z1,
        VERTEX, x0, y1, z0,
        
        END
    ]
         
    cmd.load_cgo(box,"boxxx")
