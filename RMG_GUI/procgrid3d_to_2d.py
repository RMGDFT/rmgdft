#  input:  pex, pey, pez are the numbers of processors for 3D real space
#           grids
#  output: nprow, npcol is the 2D processor grids for matrix operations
#            nprow * npcol = pex *pey * pez


#from Prime import *

def procgrid3d_to_2d(pex, pey, pez):
    
    p = Prime()
    p1 = p.factorize(pex)
    p2 = p.factorize(pey)
    p3 = p.factorize(pez)
    ptot = sorted(p1 + p2 + p3, reverse=True)

    
    nprow = 1
    npcol = 1
    for i in range(len(ptot)/2):
        npcol *= ptot[2*i]
        nprow *= ptot[2*i+1]

    if (len(ptot) - len(ptot)/2 * 2 ==1 ):
        npcol *= ptot[len(ptot)-1]

    return nprow, npcol
