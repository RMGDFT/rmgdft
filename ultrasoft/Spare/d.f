c************************** SVN Revision Information **************************
c **    $Id$    **
c******************************************************************************
 


	integer pack_doublesf
	integer idx, b(300)
	real*8 a(100), c(100)
	a(3) = 17.1246d0
	idx = pack_doublesf(a, b, 10, 5)
	call unpack_doublesf(c, b, 10, 5)
	write(6,100)c(3)
 100    format(2x,d18.12)
	end

